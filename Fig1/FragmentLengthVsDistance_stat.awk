{
	if($15>=0 && $15<=180){
		if($12=="+"){
			if($9<=$3)
				print $5";"$15
			else
				print $5";"$15*(-1)
			}
		if($12=="-"){
			if($9<=$3)
				print $5";"$15*(-1)
			else
				print $5";"$15
			}
		}
}