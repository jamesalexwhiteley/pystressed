# Author: James Whiteley (github.com/jamesalexwhiteley)

class SectionForce():
    """
    Forces and sections properties.

    """
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)