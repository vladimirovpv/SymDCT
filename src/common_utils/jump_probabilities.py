"""Manage jump probability symbols."""

from __future__ import annotations

from sympy import symbols


class JumpProbabilities:
    """Manage symbolic jump probability variables.

    Maintains a dictionary of symbolic variables (W0, W1, ...) representing
    jump probabilities, with support for customization (set to zero,  equate pairs).
    """

    def __init__(self, number_of_jumps: int):
        """Initialize jump probabilities dictionary.

        Parameters
        ----------
        number_of_jumps
            Total number of distinct jumps in the network.
        """
        self.w = {
            i: symbols(f"W{i}", real=True)
            for i in range(number_of_jumps)
        }

    def customize_jump_probabilities(self, customizations: dict) -> None:
        """Apply customizations to jump probabilities.

        Parameters
        ----------
        customizations
            Dict mapping jump indices to customizations:
            - index: other_index → set w[index] = w[other_index] (equate two jumps)
            - index: None → set w[index] = 0 (disable jump)

        Example
        -------
        >>> jp = JumpProbabilities(5)
        >>> jp.customize_jump_probabilities({2: 1, 3: None})
        Sets w[2] = w[1] and w[3] = 0.
        """
        for idx, target in customizations.items():
            if target is not None:
                # Equate to another jump's probability
                self.w[idx] = self.w[target]
            else:
                # Set to zero
                self.w[idx] = 0

    def get_jump_probabilities(self) -> dict:
        """Return the jump probabilities dictionary."""
        return self.w

