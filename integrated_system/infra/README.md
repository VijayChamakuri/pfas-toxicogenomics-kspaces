# AWS deployment design

This directory defines an approval-gated AWS target for the integrated PFAS analysis system. It has not been deployed, benchmarked in AWS, or granted cloud credentials. Local scientific results do not depend on AWS.

## Architecture

The source diagram is [architecture.mmd](architecture.mmd).

Amazon S3 stores versioned, access-controlled input snapshots and derived artifacts. Amazon ECR stores images identified by immutable commit tags and digests. ECS Fargate runs the stateless analysis-planning API behind an existing HTTPS Application Load Balancer. CloudWatch receives operational logs, metrics, and alarms. Secrets Manager is reserved for provider credentials when an approved integration needs them. A bounded SageMaker job is optional for managed model adaptation or batch scientific execution when local resources are insufficient.

The current service has no justified transactional persistence requirement, so the design does not provision RDS. It also avoids Kubernetes and Step Functions because the present system has one stateless service and explicitly launched batch jobs rather than a large persistent workflow graph.

## Data and release flow

1. Fast CI validates formatting, linting, types, unit and integration tests, fixture notebook execution, schema checks, workflow grounding, abstention, model smoke tests, and container builds.
2. Expensive full-data analysis and model adaptation run through separately approved jobs. Their outputs include data versions, image digests, configuration, seeds, metrics, and review status.
3. An approval-gated release job uses GitHub OIDC to push an image tagged with the commit SHA to ECR and register an ECS task-definition revision.
4. An existing HTTPS load balancer routes requests to healthy Fargate tasks in private subnets.
5. A pinned AWS CLI initialization container copies one approved, versioned S3 prefix to an ephemeral task volume. It exits before the API starts. The API mounts that volume read-only.
6. The service validates each request, retrieves evidence, emits a typed plan, and requires the declared review checkpoint before allowlisted execution.
7. Results and evaluation records are stored as versioned artifacts. Operational health and scientific validity remain separate release gates.

## Local validation

```bash
cd pfas-toxicogenomics-kspaces
docker build --target api -f integrated_system/docker/Dockerfile -t pfas-workflow-api:local .
docker run --rm -p 8000:8000 -v "$PWD/..:/data:ro" \
  pfas-workflow-api:local
curl --fail http://localhost:8000/health

cd integrated_system/infra/terraform
terraform init -backend=false
terraform fmt -check -recursive
terraform validate
```

Copy `terraform.tfvars.example` to the ignored `terraform.tfvars` only when preparing a reviewed plan. Replace every placeholder with approved infrastructure identifiers, then run:

```bash
terraform plan -var-file=terraform.tfvars
```

`terraform plan` does not create infrastructure. Do not run `terraform apply` without explicit approval, a cost review, and a security review.

## Security boundary

- Fargate tasks receive no public IP and accept port 8000 only from the existing load-balancer security group.
- Containers run as UID 10001 with a read-only root filesystem and contain no credentials.
- The task role has read-only access to one artifact bucket. Add write access only to an explicit versioned result prefix after the application supports cloud writes.
- GitHub deployment uses a branch-scoped or environment-scoped OIDC role, never long-lived AWS access keys.
- Provider tokens belong in Secrets Manager and are injected only into the process that requires them.
- S3 blocks public access, encrypts objects at rest, and enables versioning.
- Terraform state may contain sensitive identifiers. Production use requires an encrypted remote backend, locking, and restricted access.
- The API currently has no authentication layer. Add an identity-aware gateway or load-balancer OIDC before nonlocal use.

For regulated data, add a customer-managed KMS key, VPC endpoints, CloudTrail data events, approved retention policies, and a formal data-classification review.

## Monitoring and quality gates

Operational monitoring includes running task count, load-balancer 5xx rate, response time, unhealthy targets, CPU, memory, restarts, and structured application errors. Alerts need an explicit destination before deployment.

Scientific and model quality are separate from uptime. Each release should record dataset object versions, code SHA, image digest, model and embedding identifiers, evidence-corpus version, prompt version, latency, abstention rate, plan validity, citation correctness, unsupported-claim rate, and subgroup error slices. Drift should trigger review, not automatic retraining.

## Cost controls

Recurring costs primarily come from the load balancer, Fargate task time, NAT traffic, logs, and optional SageMaker jobs. Before deployment:

1. estimate region-specific cost with the AWS Pricing Calculator;
2. configure AWS Budgets and ownership tags;
3. bound CloudWatch retention and ECR image count;
4. set nonproduction service count to zero when idle;
5. define maximum runtime and automatic shutdown for every managed training job;
6. decide whether VPC endpoints reduce sustained NAT cost.

## Rollback and incident evidence

ECS uses a deployment circuit breaker with automatic rollback. A manual rollback selects the last healthy task-definition revision whose ECR digest and evaluation artifacts are known, then waits for service stability. Preserve the failed image, request logs, input versions, and model artifacts until the incident is understood. Restore data through S3 versioning. Revert infrastructure through source control and a reviewed Terraform plan, not ad hoc console edits.

## Terraform scope and prerequisites

The configuration creates an artifact bucket, ECR repository, ECS cluster and service, task roles, log group, security group, and a health alarm. It requires existing private subnets, an HTTPS load-balancer target group, and the load balancer's security group. Production also requires an encrypted remote state backend, OIDC deployment role, DNS, TLS certificate, authentication layer, alert destination, and approved network-egress policy. These organization-specific resources are not guessed by the starter configuration.
