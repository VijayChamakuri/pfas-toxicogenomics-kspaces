variable "aws_region" {
  type    = string
  default = "us-east-1"
}
variable "project_name" {
  type    = string
  default = "pfas-workflow"
}
variable "environment" {
  type    = string
  default = "dev"
}
variable "vpc_id" {
  type        = string
  description = "Existing governed VPC ID"
}
variable "private_subnet_ids" {
  type        = list(string)
  description = "At least two private subnet IDs"
}
variable "alb_security_group_id" {
  type        = string
  description = "Existing ALB security group ID"
}
variable "target_group_arn" {
  type        = string
  description = "Existing HTTPS ALB target group ARN"
}
variable "image_tag" {
  type        = string
  description = "Immutable image tag, preferably a commit SHA"
}
variable "desired_count" {
  type    = number
  default = 1
}
