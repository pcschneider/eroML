from eroML.ensemble import multi_fits_support
from eroML.tile import file4
import numpy as np

@multi_fits_support(2)
def props2asso(e, asso={}, props=None):
    """
    Properties for associations
    
    Example::
    
      x = props2asso(fn,asso={"srcID":["ML24537","ML27261"], "srcID_NN":["3843400348669240192","3843576060075670016"]})#, props="Fx")
    
    Parameters
    ----------
    e : Ensemble or filename
    asso : Two(!)-key dictionary
      The two keys are the column names for the associations. The value contain arrays/lists with the IDs
    props : array / list
      Properties to return. If None, return all properties
    
    Returns
    -------
    properties : dict
    
    """
    id_col0, id_col1 = list(asso.keys())[0],list( asso.keys())[1]
    arr = e.to_array((id_col0, id_col1, "srcID"), array_type="array")
    ens = [str(a)+"&"+str(b) for a, b in zip(arr[0], arr[1])]
    ser = [str(a)+"&"+str(b) for a, b in zip(asso[id_col0], asso[id_col1])]
    ii = np.where(np.in1d(ens, ser))[0]
    srcIDs = arr[2][ii]    
    e.keep(srcIDs)
    if props == None:
        props = e.known_cols
    ret = e.to_array(props, array_type = "dict")
    return ret

if __name__ == "__main__":
    fn = file4("major", cconfig="eFEDS_EDR3.ini")
    print("Using ",fn)
    x =props2asso(fn,asso={"srcID":["ML24537","ML27261"], "srcID_NN":["3843400348669240192","3843576060075670016"]})#, props="Fx")
    print(x)
    print(x["Fx"])
