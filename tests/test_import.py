"""
This is a general test for importing the tool for now
"""

import chemical_perception

try:
    import openeye
    from openeye import oechem
    TEST_OE = True
except:
    print("Tests will not be performed with openeye tools")
    TEST_OE = False

def test_e_test():
    print("This is running a function inside test_imports.py")
    if TEST_OE:
        print("and openeye was imported successfully")
