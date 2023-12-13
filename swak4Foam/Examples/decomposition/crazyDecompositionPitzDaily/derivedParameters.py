if numberOfProcessors < 2:
    from PyFoam.Error import error
    error("This case only makes sense with more than one processor")
