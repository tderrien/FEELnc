use ExtUtils::MakeMaker;

WriteMakefile(
        'NAME' => 'FEELnc',
        'VERSION_FROM' => 'lib/FEELnc.pm',
        'PREREQ_PM' =>{Bio::Tools::GFF => 0, Parallel::ForkManager => 1.06, Bio::SeqFeature::Generic=> 0},
        'AUTHOR' =>     ' Fabrice Legeai ; Audrey David ; Thomas Derrien'
);
