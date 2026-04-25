"""Shared pytest fixtures for ScanITD tests."""

from scanitd.base.basic import Event, MicroRegion


def make_event() -> Event:
    """Create a representative event used by writer/utility tests."""
    return Event.new(
        event_type="TDUP",
        event_id=("chr1", 99, 12, "ACTGACTGACTG", MicroRegion("+AA")),
        oao=8,
        ao=5,
        dp=20,
        ref_allele="A",
        alt_allele="<TDUP>",
    )
