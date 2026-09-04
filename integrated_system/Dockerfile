FROM python:3.13-slim
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .
ENV PFAS_PROJECT_ROOT=/data
EXPOSE 8000
CMD ["uvicorn", "pfas_workflow.api:app", "--host", "0.0.0.0", "--port", "8000"]
