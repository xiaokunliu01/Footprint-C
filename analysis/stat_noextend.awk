{
	n=n+1
	if(($14!=-1 && $14<=0) || ($28<=0 && $28!=-1)){
		sig=sig+1
		}
	if((($14==-1 || $14>0) && ($28<=0 && $28!=-1)) || (($14!=-1 && $14<=0) && ($28>0 || $28==-1))){
		onlysig=onlysig+1
		}
	if(($14!=-1 && $14<=0) && ($28<=0 && $28!=-1)){
		both=both+1
		if($1==$15){
			both_cis=both_cis+1
			if($21-$7>1000)
				both_cis_1k=both_cis_1k+1
			}
		}
	if(($1==$15) && ($14!=-1 && $14<=0) && ($28<=0 && $28!=-1) && ($21-$7>1000)){
		if($11=="+" && $25=="-"){
			con=con+1
			if($5=="+" && $19=="-"){
				con_FR=con_FR+1
				}
			if($5=="+" && $19=="+"){
				con_FF=con_FF+1
				}
			if($5=="-" && $19=="+"){
                        	con_RF=con_RF+1
                        	}
			if($5=="-" && $19=="-"){
                        	con_RR=con_RR+1
                        	}
			}
		if($11=="+" && $25=="+"){
                	tandemF=tandemF+1
			if($5=="+" && $19=="-"){
                        	tandemF_FR=tandemF_FR+1
                        	}
                	if($5=="+" && $19=="+"){
                        	tandemF_FF=tandemF_FF+1
                        	}
                	if($5=="-" && $19=="+"){
                        	tandemF_RF=tandemF_RF+1
                        	}
                	if($5=="-" && $19=="-"){
                        	tandemF_RR=tandemF_RR+1
                        	}
                	}
		if($11=="-" && $25=="+"){
                	div=div+1
			if($5=="+" && $19=="-"){
                        	div_FR=div_FR+1
                        	}
                	if($5=="+" && $19=="+"){
                        	div_FF=div_FF+1
                        	}
                	if($5=="-" && $19=="+"){
                        	div_RF=div_RF+1
                        	}
                	if($5=="-" && $19=="-"){
                        	div_RR=div_RR+1
                        	}
                	}
		if($11=="-" && $25=="-"){
                	tandemR=tandemR+1
			if($5=="+" && $19=="-"){
                        	tandemR_FR=tandemR_FR+1
                        	}
                	if($5=="+" && $19=="+"){
                        	tandemR_FF=tandemR_FF+1
                        	}
                	if($5=="-" && $19=="+"){
                        	tandemR_RF=tandemR_RF+1
                        	}
                	if($5=="-" && $19=="-"){
                        	tandemR_RR=tandemR_RR+1
                        	}
                	}
		}
}
END{print("pairs\t"n"\n""single\t"sig"\n""only single\t"onlysig"\n""both\t"both"\n""both cis\t"both_cis"\n""both cis (> 1k)\t"both_cis_1k"\n\n"" \tCon\tTandemF\tDiv\tTandemR\n""total\t"con"\t"tandemF"\t"div"\t"tandemR"\n""FR\t"con_FR"\t"tandemF_FR"\t"div_FR"\t"tandemR_FR"\n""FF\t"con_FF"\t"tandemF_FF"\t"div_FF"\t"tandemR_FF"\n""RF\t"con_RF"\t"tandemF_RF"\t"div_RF"\t"tandemR_RF"\n""RR\t"con_RR"\t"tandemF_RR"\t"div_RR"\t"tandemR_RR"\n")}
