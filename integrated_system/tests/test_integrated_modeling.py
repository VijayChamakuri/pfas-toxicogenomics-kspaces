from pathlib import Path

from pfas_workflow.modeling import ModelingContract, load_router_metadata, router_disclosure

ROOT = Path(__file__).parents[1]


def test_shared_modeling_contract_is_backend_neutral_and_scoped():
    contract = ModelingContract()
    assert contract.supported_backends == ("pytorch", "tensorflow")
    assert contract.seed == 2026
    assert "implementation verification" in contract.scientific_scope
    assert "probabilities" in contract.uncertainty_outputs


def test_adapted_router_metrics_are_reported_without_quality_claim():
    path = ROOT / "llm_system" / "training_output" / "metrics.json"
    metadata = load_router_metadata(path)
    assert metadata.model_name == "HuggingFaceTB/SmolLM2-135M-Instruct"
    assert metadata.training_steps == 8
    assert metadata.held_out_loss > 0
    assert "not approved" in metadata.evaluation_status
    assert router_disclosure(path)["selected_router"] == "deterministic evidence rules"
