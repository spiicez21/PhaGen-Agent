# PhaGen Kubernetes Deployment

This directory contains Kubernetes manifests for deploying PhaGen to a production cluster (AWS EKS, Azure AKS, GCP GKE, or self-hosted K8s).

## Architecture

- **API Pods**: FastAPI backend with health checks, auto-scaling (2-10 replicas)
- **Celery Workers**: Distributed job processing (1-8 replicas)
- **Celery Beat**: Periodic task scheduler (1 replica)
- **Redis**: Cache and Celery broker/backend
- **Ingress**: NGINX ingress with TLS termination

## Prerequisites

1. **Kubernetes cluster** (v1.24+)
2. **kubectl** configured to access your cluster
3. **Helm** (optional, for nginx-ingress and cert-manager)
4. **Container registry** with phagen/api:latest image

## Quick Start

### 1. Create namespace and secrets

```bash
kubectl apply -f namespace.yaml

# Edit configmap.yaml and secrets
kubectl apply -f configmap.yaml
```

**Important**: Update `phagen-secrets` in `configmap.yaml` with your actual credentials:
- `DATABASE_URL`: Connection string for PostgreSQL (Neon/RDS/CloudSQL)
- `S3_ACCESS_KEY` / `S3_SECRET_KEY`: AWS S3 or MinIO credentials
- `S3_ENDPOINT_URL`: S3 endpoint (leave blank for AWS, set for MinIO)

### 2. Deploy services

```bash
# Deploy Redis first
kubectl apply -f redis-deployment.yaml

# Deploy API
kubectl apply -f api-deployment.yaml

# Deploy Celery workers
kubectl apply -f celery-deployment.yaml

# Enable autoscaling
kubectl apply -f hpa.yaml

# Configure ingress (requires nginx-ingress-controller)
kubectl apply -f ingress.yaml
```

### 3. Verify deployment

```bash
# Check pod status
kubectl get pods -n phagen

# Check services
kubectl get svc -n phagen

# View API logs
kubectl logs -n phagen deployment/phagen-api -f

# Check autoscaling status
kubectl get hpa -n phagen
```

## Scaling

### Manual scaling

```bash
# Scale API pods
kubectl scale deployment phagen-api -n phagen --replicas=5

# Scale Celery workers
kubectl scale deployment phagen-celery-worker -n phagen --replicas=4
```

### Automatic scaling

HPA (Horizontal Pod Autoscaler) is configured to auto-scale based on:
- **API**: 2-10 replicas (CPU >70%, Memory >80%)
- **Celery Workers**: 1-8 replicas (CPU >75%, Memory >85%)

## Ingress Setup

### Install nginx-ingress-controller

```bash
helm repo add ingress-nginx https://kubernetes.github.io/ingress-nginx
helm install nginx-ingress ingress-nginx/ingress-nginx -n ingress-nginx --create-namespace
```

### Install cert-manager for TLS

```bash
helm repo add jetstack https://charts.jetstack.io
helm install cert-manager jetstack/cert-manager --namespace cert-manager --create-namespace --set installCRDs=true

# Create Let's Encrypt issuer
kubectl apply -f - <<EOF
apiVersion: cert-manager.io/v1
kind: ClusterIssuer
metadata:
  name: letsencrypt-prod
spec:
  acme:
    server: https://acme-v02.api.letsencrypt.org/directory
    email: admin@example.com
    privateKeySecretRef:
      name: letsencrypt-prod
    solvers:
    - http01:
        ingress:
          class: nginx
EOF
```

Update `ingress.yaml` with your domain and apply.

## Monitoring

### View metrics

```bash
# API metrics (Prometheus format)
kubectl port-forward -n phagen svc/phagen-api 8000:8000
curl http://localhost:8000/metrics

# Check health
curl http://localhost:8000/health
curl http://localhost:8000/ready
```

### Logs

```bash
# API logs
kubectl logs -n phagen -l component=api --tail=100 -f

# Celery worker logs
kubectl logs -n phagen -l component=celery-worker --tail=100 -f

# Redis logs
kubectl logs -n phagen -l component=redis --tail=100 -f
```

## Production Checklist

- [ ] Update secrets in `configmap.yaml` with real credentials
- [ ] Configure external PostgreSQL (Neon/RDS/CloudSQL) in `DATABASE_URL`
- [ ] Set up S3 buckets and configure credentials
- [ ] Update ingress host to your domain
- [ ] Configure TLS certificates (Let's Encrypt or custom)
- [ ] Set resource limits based on load testing
- [ ] Enable pod disruption budgets for high availability
- [ ] Configure backup strategy for Redis and Postgres
- [ ] Set up monitoring (Prometheus + Grafana)
- [ ] Configure log aggregation (ELK/Loki)
- [ ] Implement network policies for security

## Resource Requirements

### Minimum (development)
- 2 vCPU, 4GB RAM total
- API: 2 pods × 250m CPU, 512Mi RAM
- Celery: 1 worker × 500m CPU, 1Gi RAM
- Redis: 100m CPU, 256Mi RAM

### Recommended (production)
- 8 vCPU, 16GB RAM total
- API: 3-5 pods × 500m CPU, 1Gi RAM
- Celery: 2-4 workers × 1 CPU, 2Gi RAM
- Redis: 500m CPU, 512Mi RAM

## Troubleshooting

### Pods not starting

```bash
kubectl describe pod -n phagen <pod-name>
kubectl logs -n phagen <pod-name>
```

### Database connection issues

```bash
# Test from inside a pod
kubectl exec -it -n phagen deployment/phagen-api -- python -c "from backend.app.database import SessionLocal; db = SessionLocal(); print('DB OK')"
```

### Redis connection issues

```bash
# Test Redis connectivity
kubectl exec -it -n phagen deployment/phagen-redis -- redis-cli ping
```

### Celery workers not processing jobs

```bash
# Check Celery worker logs
kubectl logs -n phagen -l component=celery-worker

# Verify Redis connectivity from worker
kubectl exec -it -n phagen deployment/phagen-celery-worker -- python -c "import redis; r = redis.from_url('redis://phagen-redis:6379/1'); print(r.ping())"
```

## Cleanup

```bash
# Delete all resources
kubectl delete namespace phagen

# Or delete individually
kubectl delete -f ingress.yaml
kubectl delete -f hpa.yaml
kubectl delete -f celery-deployment.yaml
kubectl delete -f api-deployment.yaml
kubectl delete -f redis-deployment.yaml
kubectl delete -f configmap.yaml
kubectl delete -f namespace.yaml
```
