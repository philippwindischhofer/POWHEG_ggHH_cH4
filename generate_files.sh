#! /bin/bash

if [[ $1 == 'btilde' ]]
then
    
    cat ../btilde.f | perl -pe 's/         call allborn/         www=www*BI_rescaling()\n         wwwtot=wwwtot*BI_rescaling()\n         call allborn/ ; s/      if\(ini\) then/      real * 8 BI_rescaling\n      external BI_rescaling\n      if(ini) then/' > btilde_gghh.f
    
    replacements_btilde=$(grep -c 'BI_rescaling' btilde_gghh.f)
    
    if [[ $replacements_btilde == '4' ]]
    then
	echo 'File btilde_gghh.f created correctly!'
    else
	echo 'There was a problem creating btilde_gghh.f, please contact the developers!!'
	rm btilde_gghh.f
    fi

elif [[ $1 == 'sigremnants' ]]
then
    
    cat ../sigremnants.f | perl -pe 's/logical pwhg_isfinite/real * 8 BI_rescaling\n      external BI_rescaling\n      logical pwhg_isfinite/ ; s/xjac=jac_over_csi\*kn_csi\*kn_csimax\*kn_jacborn\*ww\*hc2/xjac=jac_over_csi\*kn_csi\*kn_csimax\*kn_jacborn\*ww\*hc2\*\n     # BI_rescaling\(\)/ ; s/     #                \*kn_jacborn\*ww\*hc2/     #                \*kn_jacborn\*ww\*hc2\*BI_rescaling\(\)/' > sigremnants_gghh.f

    replacements_sigremnants=$(grep -c 'BI_rescaling' sigremnants_gghh.f)
    
    if [[ $replacements_sigremnants == '4' ]]
    then
	echo 'File sigremnants_gghh.f created correctly!'
    else
	echo 'There was a problem creating sigrenmants_gghh.f, please contact the developers!!'
	#rm sigremnants_gghh.f
    fi

else
    echo 'File name not recognized, abort.'
    exit
fi
