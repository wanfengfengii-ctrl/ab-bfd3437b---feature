FROM python:3.11-slim

WORKDIR /srv

# 服务仅使用 Python 标准库，无第三方依赖
COPY app/ ./app/
COPY tests/ ./tests/
COPY verify.py ./verify.py

ENV API_PORT=8080 \
    PYTHONUNBUFFERED=1

EXPOSE 8080

# 默认运行 API；verify 服务在 compose 中以 `python verify.py` 覆盖启动命令
CMD ["python", "app/main.py"]
