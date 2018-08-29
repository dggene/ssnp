Channel
    .from( 'alpha,beta,gamma\n10,20,30\n70,80,90' )
    .splitText(file:true)
    .into{rows0;rows1}




process first_line{
	input:
		val header from rows0.first()
		val row from rows1
	
	script:
	"""
		echo $header
		echo $row
	"""
}