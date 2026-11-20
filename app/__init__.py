"""Flask 应用工厂与路由。"""

from __future__ import annotations

import os

from flask import (
    Flask,
    abort,
    flash,
    redirect,
    render_template,
    request,
    url_for,
)

from . import services
from .db import init_app
from .normalization import ValidationError


def create_app() -> Flask:
    app = Flask(__name__)
    app.config["DATABASE"] = os.environ.get(
        "DATABASE_PATH", "/data/variants.db"
    )
    app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY", "dev-secret-change-me")
    app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # 16 MB

    init_app(app)

    @app.get("/")
    def index():
        data = services.get_dashboard()
        return render_template("index.html", **data)

    @app.post("/reference")
    def upload_reference():
        file = request.files.get("reference")
        if file is None or not file.filename:
            flash("请选择参考序列文件。", "error")
            return redirect(url_for("index"))
        try:
            result = services.upload_reference(file.read())
        except ValidationError as exc:
            flash(f"参考上传失败，当前结果未改变：{exc}", "error")
            return redirect(url_for("index"))
        flash(
            f"参考已更新（长度 {result['length']}）。重新校验 "
            f"{len(result['revalidated'])} 个样本通过，"
            f"{len(result['failed'])} 个样本标记为待重传。",
            "success",
        )
        for row in result["failed"]:
            flash(f"样本「{row['name']}」待重传：{row['error']}", "warning")
        return redirect(url_for("index"))

    @app.post("/samples")
    def upload_sample():
        name = request.form.get("name", "")
        file = request.files.get("sample")
        if file is None or not file.filename:
            flash("请选择样本 TSV 文件。", "error")
            return redirect(url_for("index"))
        try:
            result = services.upload_sample(name, file.read())
        except ValidationError as exc:
            flash(f"样本上传失败，当前结果未改变：{exc}", "error")
            return redirect(url_for("index"))
        action = "替换旧版本" if result["replaced"] else "新建"
        flash(
            f"样本「{result['sample']}」{action}成功，"
            f"归一化后 {result['count']} 条变异。",
            "success",
        )
        return redirect(url_for("index"))

    @app.post("/samples/<int:sample_id>/delete")
    def delete_sample(sample_id: int):
        services.delete_sample(sample_id)
        flash("样本已删除，相关变异计数已更新。", "success")
        return redirect(url_for("index"))

    @app.get("/variants/<int:pos>/<path:alleles>")
    def variant_detail(pos: int, alleles: str):
        if alleles.count(">") != 1:
            abort(404)
        ref, alt = alleles.split(">", 1)
        data = services.get_variant_detail(pos, ref, alt)
        if data is None:
            abort(404)
        return render_template("variant.html", variant=data)

    @app.get("/healthz")
    def healthz():
        return {"status": "ok"}

    @app.errorhandler(404)
    def not_found(_err):
        return render_template("404.html"), 404

    @app.errorhandler(413)
    def too_large(_err):
        flash("上传文件超过 16 MB 限制。", "error")
        return redirect(url_for("index"))

    return app
