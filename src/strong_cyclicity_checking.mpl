strong_cyclicity_checking := module()
    option package;
    local prolong, prolong_matrix_one, show;
    export is_strong_cyclic, CV_is_cyclic, CV_dv, CV_companion_matrix, MAX_RECURSION_DEPTH, VERBOSE;

$include "cyclic.mm"
$include "prolong.mm"
$include "strong_cyclic.mm"

end module;
