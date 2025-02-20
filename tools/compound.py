"""
This file describes a compound class
"""

import re
from collections import defaultdict

class Compound:
    """
    A class to describe compounds in a metabolic network.
    """
    def __init__(self, items, **kwargs):
        self.r_id          = items[0]          # An id string
        self.name          = items[1]
        self.compartment   = items[2]
        self.concentration = items[3]
        self.charge        = items[4]
        self.formula       = items[5]
        self.dict          = {}
        for k, v in kwargs.items():
            self.dict[k] = v

    @classmethod
    def from_text(cls, line):
        """
        Create a compound based on a csv line containing the information.
        """
        items = line.split(',')
        items = items[:2] + [items[0][1]] + items[2:]
        return cls(items)

    def atoms(self):
        """
        Return the formula as a dictionary of 'Element', ammount pairs
        """
        pattern = r"([A-Z][a-z]?)(\d*(?:\.\d+)?)"

        element_counts = defaultdict(float)

        for element, count in re.findall(pattern, self.formula):
            count = float(count) if count else 1.0  # Convert to float, default to 1.0 if missing
            element_counts[element] += count  # Sum up duplicate elements

        return dict(element_counts)
