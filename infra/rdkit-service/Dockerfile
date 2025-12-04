FROM mambaorg/micromamba:1.5.8

RUN apt-get update && apt-get install -y --no-install-recommends curl && rm -rf /var/lib/apt/lists/*

WORKDIR /service

COPY environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml && \
    micromamba clean --all --yes
ENV PATH=/opt/conda/envs/rdkit-service/bin:$PATH \
    CONDA_DEFAULT_ENV=rdkit-service

COPY app.py /service/app.py

EXPOSE 8081

CMD ["/bin/bash", "-lc", "micromamba run -n rdkit-service uvicorn app:app --host 0.0.0.0 --port 8081"]
