StrongCyclicityChecking := module()
    option package;
    local prolong, remove_o_notation, max_degree, show, find_roots, is_strong_cyclic_process;
    export is_strong_cyclic, CV_is_cyclic, CV_dv, CV_companion_matrix;
    global RESULT_MATRIX, VERBOSE, MAX_RECURSION_DEPTH;

$include "cyclic.mm"
$include "polynom.mm"
$include "strong_cyclic.mm"

end module;
