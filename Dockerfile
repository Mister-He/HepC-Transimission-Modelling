# =============================================================================
# Dockerfile — NumPyro/JAX 推断环境
#   docker build -t hepc-numpyro:latest .
#   docker run --rm hepc-numpyro:latest                # 运行全部测试
#   docker run --rm -v "$PWD":/workspace hepc-numpyro:latest \
#       python scripts/run_nuts.py --num-warmup 500 --num-samples 1000
# =============================================================================
FROM python:3.13-slim

WORKDIR /workspace

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV MPLCONFIGDIR=/tmp/mplcache

CMD ["python", "-m", "pytest", "tests/unit", "-v"]
