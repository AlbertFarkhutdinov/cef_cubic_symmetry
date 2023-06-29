"""The module contains CEF class."""


from pretty_repr import RepresentableObject

from core.sample import Sample


class BaseInteraction(RepresentableObject):
    """
    A base class for interactions of a sample with fields.

    Parameters
    ----------
    sample : core.Sample
        A sample interacting with fields.

    """

    def __init__(self, sample: Sample) -> None:
        """Initialize self. See help(type(self)) for accurate signature."""
        self.sample = sample
