"""Tests for VCF writer helpers."""

from scanitd.base.basic import Event, MicroRegion
from scanitd.writer.vcf_writer import get_vcf_features_from_event, obtain_reference_from_bam_header, vcf_feature_transformer


def test_obtain_reference_from_bam_header_finds_known_aligner() -> None:
    bam_header = {"PG": [{"ID": "bwa", "CL": "bwa mem ref.fa in.bam"}]}

    assert obtain_reference_from_bam_header(bam_header) == "bwa mem ref.fa in.bam"


def test_vcf_feature_transformer_has_gt_field() -> None:
    event = Event.new(
        event_type="TDUP",
        event_id=("chr1", 99, 12, "ACTGACTGACTG", MicroRegion("+AA")),
        oao=8,
        ao=5,
        dp=20,
        ref_allele="A",
        alt_allele="<TDUP>",
    )

    feature_dict = get_vcf_features_from_event(event)
    transformed = vcf_feature_transformer(feature_dict, event_id=7)

    assert transformed[0] == "chr1"
    assert transformed[2] == "7"
    assert transformed[8] == "GT"
    assert transformed[9] == "0/1"
    assert "SVTYPE=TDUP" in transformed[7]
    assert "INSSEQ=AA" in transformed[7]
