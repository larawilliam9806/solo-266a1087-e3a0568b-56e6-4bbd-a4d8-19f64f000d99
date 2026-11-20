FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DATABASE_PATH=/data/variants.db

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app
COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh && useradd -r -u 1000 appuser && mkdir -p /data && chown -R appuser /data

VOLUME ["/data"]
EXPOSE 8000
USER appuser

ENTRYPOINT ["/entrypoint.sh"]
