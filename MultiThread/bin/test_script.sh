	
	sampleData[0]=["/mirror/Test/TestData/TestSamples/sample1000_3.txt"]
	sampleData[1]=["/mirror/Test/TestData/TestSamples/sample10000_3.txt"]
	sampleData[2]=["/mirror/Test/TestData/TestSamples/sample100000_3.txt"]
	sampleData[3]=["/mirror/Test/TestData/TestSamples/sample1000000_3.txt"]
	sampleData[4]=["/mirror/Test/TestData/TestSamples/sample1000_4.txt"]
	sampleData[5]=["/mirror/Test/TestData/TestSamples/sample10000_4.txt"]
	sampleData[6]=["/mirror/Test/TestData/TestSamples/sample100000_4.txt"]
	sampleData[7]=["/mirror/Test/TestData/TestSamples/sample1000000_4.txt"]
	sampleData[8]=["/mirror/Test/TestData/TestSamples/sample1000_5.txt"]
	sampleData[9]=["/mirror/Test/TestData/TestSamples/sample10000_5.txt"]
	sampleData[10]=["/mirror/Test/TestData/TestSamples/sample10000_5.txt"]
	sampleData[11]=["/mirror/Test/TestData/TestSamples/sample100000_5.txt"]
	sampleData[12]=["/mirror/Test/TestData/TestSamples/sample1000_6.txt"]
	sampleData[13]=["/mirror/Test/TestData/TestSamples/sample10000_6.txt"]
	sampleData[14]=["/mirror/Test/TestData/TestSamples/sample100000_6.txt"]
	sampleData[15]=["/mirror/Test/TestData/TestSamples/sample1000000_6.txt"]
	finalResult[0]=0
	dataSet[0]=[0]
	row=3
	cols=1000
        max_threads=16
	data_source="/mirror/Test/TestData/TestSamples/sample1000_3.txt"
        for j in {1..$max_threads}
	do		
		for k in {1..3}
		do
			./hFastICA -C cubic -I 500 -c $row -g identity -i /mirror/Test/TestData/TestSamples/sample1000_3.txt -r 200 -s $cols -t 0.0000001 -if text -C cubic -ow ./result1.txt -th $max_threads
		done
	done
	count=0
	step=0
	
	IFS=$'\n' read -d '' -r -a lines < testResult.txt
	for (( c=0; c<${#lines[@]}; c++ ))
	do
		dataSet[c]=${lines[c]}
	done
	
	for (( i=0; i<$((${#lines[@]}/15)); i++ ))
	do
		echo ${dataSet[c]}
		for (( j=0; j<5; j++ ))
		do
			 finalResult[count]=$(echo ${dataSet[step]} + ${dataSet[step+5]} + ${dataSet[step+10]} | bc)
			 fg=${finalResult[count]}
			 fg=$(echo " scale=6; $fg / 3" | bc -l)
			 finalResult[count]=$fg
			 (( step += 1 ))
			 (( count += 1 ))
		done
		(( step += 10 ))
	done
	echo "================================Test result=========================================="
	step=0
    for (( c=0; c<$((count/5)); c++ ))
	do
		echo ${finalResult[step]} ${finalResult[step+1]} ${finalResult[step+2]} ${finalResult[step+3]} ${finalResult[step+4]} 
		(( step += 5 ))
	done




























