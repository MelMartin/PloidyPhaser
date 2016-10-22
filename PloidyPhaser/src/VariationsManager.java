import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class VariationsManager {
	VCFParser vm;
	int varMatIndexesSize=0;
	HashMap<VariationData,Integer> varMatIndexes=new HashMap<VariationData,Integer>();//keeps track of vcf positions<key> that contains a 'ref' 
	//different than an 'alt', and their Index order <value>,in order to build the matrix,
	HashMap<String,VariationData> varDatIds=new HashMap<String,VariationData>();//keeps track of vcf identifiers<key> and the 
	//corresponding VariationData 
	List<VariationData> variations = new ArrayList<VariationData>();
	int[][]varMat;


	public VariationsManager(VCFParser vm){
		this.vm=vm;
	}


	public void registerVariation(VariationData vd){

		vd.setMatrixIndex(++varMatIndexesSize);
		varMatIndexes.put(vd,vd.matrixIndex);//stores reference	allele	
		varDatIds.put(vd.id, vd);
		//System.out.print(vd.id+"/"+vd.matrixIndex);
		VariationData newVD=vd.makeAlterativeVarDataFromRef(vd);	
		newVD.setMatrixIndex(varMatIndexesSize);
		varMatIndexes.put(newVD,newVD.matrixIndex);//stores alternative allele	
		varDatIds.put(newVD.id, newVD);
		//System.out.println("   "+newVD.id+"/"+newVD.matrixIndex);
		//System.out.print("I THINK THESE TWO HASHMAP ARE NOT CORRECT FOR WHAT I NEED... CHANGE TO GET TE INDEX, THE VARDAT AND THE ID CORRECTLY");

	}


	public void fillVariantMatrix() {

		//fill VariationsPos vector (vector containing all variations positions ordered in a vector)
		int[] variationsPos=new int[variations.size()];
		for (int v=0;v<variations.size();v++){
			variationsPos[v]=variations.get(v).pos;
		}	

		int vpInd=0;//index pointing to the variationPos vector where the read maps its first variation

		SAMFileReader inputSam;
		for (int bi=0;bi<PloidyPhaser.bamPaths.size();bi++){//for each bam file
			inputSam = new SAMFileReader(PloidyPhaser.bamPaths.get(bi));
			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter = inputSam.iterator();
			int ct = 0;

			String refName = "";// for debugging a bad line in the .sam file
			int alStart;
			int alEnd;
			String libraryName="."+PloidyPhaser.bamPaths.get(bi).getName().substring(0,PloidyPhaser.bamPaths.get(bi).getName().lastIndexOf('.'));
			String readName="";
			String readSeq="";
			int currentVpInd=0;//index of the current variation

			while (iter.hasNext() && ct<500 ) {//iterates the sam file (for each read)
				SAMRecord rec = iter.next();
				readName=rec.getReadName()+libraryName;
				refName = rec.getReferenceName();
				alStart = rec.getAlignmentStart();
				alEnd = rec.getReadLength()+alStart;
				readSeq=rec.getReadString();

				//identify  variation signature of the read
				if (alStart<10000 && ct<500 ){					
					//identify start var_position covered by read (first var):

					//chek that vpInd is correctly positioned	
					while(variationsPos[vpInd]<alStart)	vpInd++;//since bam is sorted, the index only increases
					currentVpInd=vpInd;//actualize current variant position index
					if(alEnd>variationsPos[currentVpInd]){
						System.out.println("ct:"+ct+" ReadIndex:"+PloidyPhaser.readIndexes.get(readName) +" st:"+alStart+" end:"+alEnd+" length:"+rec.getReadLength()+" cig:"+rec.getCigar()+" s:"+readSeq);
						//System.out.println( "vpInd:"+vpInd+" currentVpInd:"+currentVpInd+" variationsPos_currentVpInd["+currentVpInd+"]:"+variationsPos[currentVpInd]);

					}

					//while the read spans positions with variation
					while(alEnd>variationsPos[currentVpInd]){

						//get all vars covered by the read
						int i=0;//index over all posible variations for the current read
						VariationData curVar=variations.get(currentVpInd+i);//currentVariation
						while(curVar.pos==variationsPos[currentVpInd]){
							solveVariation(curVar, alStart,readSeq,rec.getReadLength(),rec.getCigar());
							i++;
							curVar=variations.get(currentVpInd+i);

						}
						currentVpInd++;
					}		

					ct++;
				}
			}
			iter.close();
			inputSam.close();
		}
	}

	//for each variation, this function compares the subsequence of the candidate read to the possible allele variations
	//and outputs the represented one (ref=1 alt=2 none=0)
	public void solveVariation(VariationData curVar,int alStart, String readSeq, int readLength,Cigar cigarString){
		int subSeqBeg=curVar.pos-alStart;//begin index on the subsequence from readSeq.
		int subSeqEnd=(curVar.pos-alStart)+curVar.ref.length();//begin index on the subsequence from readSeq.

		if((subSeqBeg+curVar.ref.length())<=readLength){
			System.out.print( " CurrVar:"+curVar.outString()+" subSeqBeg:"+subSeqBeg);	
			CigarCounter cigCount=new CigarCounter(cigarString.toString());
			for (int i=0;i<subSeqBeg;i++){
				cigCount.next();
			}
			subSeqBeg=cigCount.beg;
			for (int i=0;i<curVar.ref.length();i++){
				cigCount.next();
			}
			subSeqEnd=cigCount.beg;
			System.out.println( " cigCount b/e:"+subSeqBeg+"/"+subSeqEnd+" query:"+readSeq.substring(subSeqBeg,subSeqEnd));
		}



	}

	private class CigarCounter{//keeps track of the way the read is mapped in order to obtain the correctly aligned read subsequence
		ArrayList<String> cigarElems;
		int beg=0;//begin index on the subsequence from readSeq.
		int end=0;//begin index on the subsequence from readSeq.
		int currentIndCigarElem=0;//current index on cigarElems
		public CigarCounter(String cigarString){
			cigarElems=splitCIGAR(cigarString);
			System.out.println( "\tcigar:"+cigarElems);
		}

		//A function to split the CIGAR value into a list of CIGAR elements 
		// with first the type of cigar alignment, then the number of corresponding position
		//ie: 60M5D40M outputs: [M, 60, D, 5, M, 40]
		private ArrayList<String> splitCIGAR(String cigarString) {
			//One cigar component is a number of any digit, followed by a letter or =
			Pattern cigarPattern = Pattern.compile("[\\d]+[a-zA-Z|=]");
			ArrayList<String> cigarElems = new ArrayList<String>();
			Matcher matcher = cigarPattern.matcher(cigarString);
			while (matcher.find()) {
				String[] unit=matcher.group().split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");	    	
				cigarElems.add(unit[1] );
				cigarElems.add(unit[0] );	    	
			}
			return cigarElems;
		}

		public void next(){
			if (beg!=-1){//this rejects all reads which contain a bad cigar Op ('S','H'or'X')
				//Could be more effective if it considers to reject only the bad Op which concerns Variation positions
				String currentOp=cigarElems.get(currentIndCigarElem);//current cigar operation (M,D,I,...)

				switch(currentOp){
				case "M"://Match
					beg++;
					break;
				case "="://equal match
					beg++;
					break;
				case "D"://deletion: don't move seqEnd
					break;
				case "I"://insertion: don't know what to do with this yet
					break;
				default://'S','H'or'X': sequence doesn't appear on alignment (don't map)->if position concerns the variation, signal to ignore the variation(returns -1?)
					beg=-1;
					break;
				}
				if (Integer.parseInt(cigarElems.get(currentIndCigarElem+1))==0 && cigarElems.size()>(currentIndCigarElem+1)){
					currentIndCigarElem=currentIndCigarElem+2;
				}
			}
		}
	}



	public  void constructVariantMatrix() {
		varMat=new int[varMatIndexes.size()][varMatIndexes.size()];

	}

	public void printOutVariations(String apendix) throws FileNotFoundException {

		System.out.println("outputCleanFile:"+vm.vcfFile.getParent()+"\\Clean\\"+vm.vcfFile.getName().substring(0, vm.vcfFile.getName().lastIndexOf('.'))+apendix+".vcf");
		vm.outputCleanFile = vm.vcfFile.getParent()+"\\Clean\\"+vm.vcfFile.getName().substring(0, vm.vcfFile.getName().lastIndexOf('.'));

		vm.outputCleanFile +=apendix+".vcf";

		PrintStream stdout = System.out;
		PrintStream myConsole = null;

		try {
			myConsole = new PrintStream(new File(vm.outputCleanFile));
			System.setOut(myConsole);
			System.out.println(vm.header);
			for (int i = 0; i < (variations.size()-1); i++) {
				System.out.println(variations.get(i).outString());
			}
			System.out.print(variations.get(variations.size()-1).outString());//last line without endline
			myConsole.close();
		} catch (Exception e) {
			System.err.println("Error trying to write outputHumanFile");
			e.printStackTrace();
		}  finally {
			if (myConsole != null) {	        	
				myConsole.close();
				System.setOut(stdout);  
			}
		}
		System.setOut(stdout);  
	}
}
