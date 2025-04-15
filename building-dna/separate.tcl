mol load pdb dsdna_raw.pdb
set all [atomselect top all]
$all moveby [vecinvert [measure center $all weight mass]]
$all moveby "0 0 20"
foreach chain {A B} {
    [atomselect top "chain $chain"] writepdb dsdna_$chain.pdb
}
