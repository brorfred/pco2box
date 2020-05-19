
from pco2_east_west import model

def test_setup():

    md = model.BoxModel(svec=730)
    md.run()
