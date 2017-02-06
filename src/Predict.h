
void C_splitnode(SEXP node, SEXP learnsample, SEXP control);

int C_get_nodeID(SEXP subtree, SEXP newinputs,
                  double mincriterion, int numobs, int varperm);
