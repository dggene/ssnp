#!/bin/env perl
package CalcHCU;
use Bio::SeqIO;

sub calc {
    shift;
    my ($seqObjs, $config) = @_;
    #return var
    my $res;
    my ($HCU_WT, $HCU_MT, $dHCU);
    my $cfFile = $config->{cfFile};
    #temp var
    my %hash;
    my %standard_codon = (
        "TTT" => 'F',
        "TTC" => 'F',
        "TTA" => 'L',
        "TTG" => 'L',
        "TCT" => 'S',
        "TCC" => 'S',
        "TCA" => 'S',
        "TCG" => 'S',
        "TAT" => 'Y',
        "TAC" => 'Y',
        "TAA" => 'B',
        "TAG" => 'B',
        "TGT" => 'C',
        "TGC" => 'C',
        "TGA" => 'B',
        "TGG" => 'W',
        "CTT" => 'L',
        "CTC" => 'L',
        "CTA" => 'L',
        "CTG" => 'L',
        "CCT" => 'P',
        "CCC" => 'P',
        "CCA" => 'P',
        "CCG" => 'P',
        "CAT" => 'H',
        "CAC" => 'H',
        "CAA" => 'Q',
        "CAG" => 'Q',
        "CGT" => 'R',
        "CGC" => 'R',
        "CGA" => 'R',
        "CGG" => 'R',
        "ATT" => 'I',
        "ATC" => 'I',
        "ATA" => 'I',
        "ATG" => 'M',
        "ACT" => 'T',
        "ACC" => 'T',
        "ACA" => 'T',
        "ACG" => 'T',
        "AAT" => 'N',
        "AAC" => 'N',
        "AAA" => 'K',
        "AAG" => 'K',
        "AGT" => 'S',
        "AGC" => 'S',
        "AGA" => 'R',
        "AGG" => 'R',
        "GTT" => 'V',
        "GTC" => 'V',
        "GTA" => 'V',
        "GTG" => 'V',
        "GCT" => 'A',
        "GCC" => 'A',
        "GCA" => 'A',
        "GCG" => 'A',
        "GAT" => 'D',
        "GAC" => 'D',
        "GAA" => 'E',
        "GAG" => 'E',
        "GGT" => 'G',
        "GGC" => 'G',
        "GGA" => 'G',
        "GGG" => 'G'
    );
    my %codon_num = (
        "TTT" => 2,
        "TTC" => 2,
        "TTA" => 6,
        "TTG" => 6,
        "TCT" => 6,
        "TCC" => 6,
        "TCA" => 6,
        "TCG" => 6,
        "TAT" => 2,
        "TAC" => 2,
        "TAA" => 3,
        "TAG" => 3,
        "TGT" => 2,
        "TGC" => 2,
        "TGA" => 3,
        "TGG" => 1,
        "CTT" => 6,
        "CTC" => 6,
        "CTA" => 6,
        "CTG" => 6,
        "CCT" => 4,
        "CCC" => 4,
        "CCA" => 4,
        "CCG" => 4,
        "CAT" => 2,
        "CAC" => 2,
        "CAA" => 2,
        "CAG" => 2,
        "CGT" => 6,
        "CGC" => 6,
        "CGA" => 6,
        "CGG" => 6,
        "ATT" => 3,
        "ATC" => 3,
        "ATA" => 3,
        "ATG" => 1,
        "ACT" => 4,
        "ACC" => 4,
        "ACA" => 4,
        "ACG" => 4,
        "AAT" => 2,
        "AAC" => 2,
        "AAA" => 2,
        "AAG" => 2,
        "AGT" => 6,
        "AGC" => 6,
        "AGA" => 6,
        "AGG" => 6,
        "GTT" => 4,
        "GTC" => 4,
        "GTA" => 4,
        "GTG" => 4,
        "GCT" => 4,
        "GCC" => 4,
        "GCA" => 4,
        "GCG" => 4,
        "GAT" => 2,
        "GAC" => 2,
        "GAA" => 2,
        "GAG" => 2,
        "GGT" => 4,
        "GGC" => 4,
        "GGA" => 4,
        "GGG" => 4
    );
    open(CF, $cfFile) or die "cannot open $cfFile,$!\n";
    while ($line = <CF>) {
        @items = split(/\t/, $line);
        ($cd, $frequency) = ($items[0], $items[2]);
        $hash{$cd} = $frequency;
    }
    close(CF);
    foreach my $seqobj (@$seqObjs) {
        $n = 0;
        @items = split(/\|/, $seqobj->display_id);
        $type = $items[2];
        $cdsLoc = $items[7];
        my $seq = $seqobj->seq;
        @codons = $seq =~ /\w{3}/g;
        if ($cdsLoc % 3 == 0) {
            $the_codon = $codons[$cdsLoc / 3 - 1];
        }
        else {
            $the_codon = $codons[int($cdsLoc / 3)];
        }
        $n = $codon_num{$the_codon};
        $frq = $hash{$the_codon};
        $HCU = $n * $frq;
        if ($type =~ /WT/) {
            $HCU_WT = $HCU;
        }
        else {
            $HCU_MT = $HCU;
        }
    }
    if ($HCU_WT eq "NA") {
        $dHCU = "NA";
    }
    $dHCU = $HCU_MT - $HCU_WT;
    $res = {
        WT   => $HCU_WT,
        MT   => $HCU_MT,
        diff => $dHCU
    };
    return $res;
}

$wt_seq_file=$ARGV[0]
$mt_seq_file=$ARGV[1]
$cf_file=$ARGV[2]



1;