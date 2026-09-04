# Security and responsible operation

## Trust boundaries

Natural-language requests, uploaded files, retrieved text, model output, environment variables, and cloud artifacts are untrusted until validated. The planner may propose a workflow, but it cannot bypass the typed schema, operation allowlist, artifact catalog, or human-review checkpoint.

## Local operation

- Run the service with read-only access to research inputs.
- Keep credentials out of notebooks, source files, images, logs, and evaluation cases.
- Bind development services to localhost unless an authenticated gateway is configured.
- Validate artifact hashes before execution and preserve the result manifest.
- Do not allow user text to become shell, Python, SQL, or path expressions.
- Restrict file access to cataloged artifacts under the configured project root.

## Model and retrieval risks

Prompt injection inside a request or retrieved document must be treated as data, not as an instruction to alter system policy. Evidence identifiers must resolve to the versioned corpus. Generated fields, citations, genes, pathways, tools, and files must be checked against their authoritative catalogs. Unsupported or conflicting evidence triggers abstention or human review.

## Cloud operation

The AWS design uses private ECS tasks, immutable ECR image tags, versioned S3 objects, bounded CloudWatch retention, and separate task and execution roles. GitHub deployment should use short-lived OIDC credentials. Runtime secrets belong in Secrets Manager and should be exposed only to the process that needs them.

Before any deployment, add authentication and authorization, define audit-log retention, configure alert destinations, review network egress, scan container images and dependencies, validate least-privilege IAM policies, and complete data-classification and cost reviews.

## Incident response

If an artifact, credential, or scientific result is suspected to be compromised:

1. stop affected execution and preserve logs;
2. revoke exposed credentials and isolate affected resources;
3. identify the exact input hash, image digest, model version, and request set;
4. invalidate derived artifacts that depend on compromised material;
5. restore a known-good version through reviewed configuration;
6. document impact on both operational integrity and scientific conclusions.

Security issues should be reported privately to the project owner. No public vulnerability-reporting address has been designated.
