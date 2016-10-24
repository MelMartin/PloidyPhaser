import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
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

			//String refName = "";// for debugging a bad line in the .sam file
			int alStart;
			int alEnd;
			String libraryName="."+PloidyPhaser.bamPaths.get(bi).getName().substring(0,PloidyPhaser.bamPaths.get(bi).getName().lastIndexOf('.'));
			String readName="";
			String readSeq="";
			int currentVpInd=0;//index of the current variation
			CigarCounter cigCount;

			while (iter.hasNext() && ct<500 ) {//iterates the sam file (for each read)
				SAMRecord rec = iter.next();//read line
				readName=rec.getReadName()+libraryName;
				//refName = rec.getReferenceName();
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
						cigCount=new CigarCounter(rec.getCigar(),readSeq);//new cigarCount for each read Sequence readSeq

						//while the read spans positions with a variation
						while(alEnd>variationsPos[currentVpInd]){

							//get all vars covered by the read
							int i=0;//index over all posible variations for the current read
							VariationData curVar=variations.get(currentVpInd+i);//currentVariation
							while(curVar.pos==variationsPos[currentVpInd]){
								solveVariation(curVar, alStart,readSeq,rec.getReadLength(),cigCount);
								i++;
								curVar=variations.get(currentVpInd+i);
							}
							currentVpInd++;
						}		
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
	public void solveVariation(VariationData curVar,int alStart, String readSeq, int readLength,CigarCounter cigCount){
		int subSeqBeg=curVar.pos-alStart;//begin index on the subsequence from readSeq.
		int subSeqEnd=(curVar.pos-alStart)+curVar.ref.length();//begin index on the subsequence from readSeq.

		if((subSeqBeg+curVar.ref.length())<=readLength){	
			System.out.println( " CurrVar:"+curVar.outString()+" subSeq:"+cigCount.getSubseq(subSeqBeg, subSeqEnd));
		}



	}

	private class CigarCounter{//keeps track of the way the read is mapped in order to obtain the correctly aligned read subsequence
		/*
		op    Description
		M    Alignment match (can be a sequence match or mismatch
		I    Insertion to the reference
		D    Deletion from the reference
		N    Skipped region from the reference
		S    Soft clip on the read (clipped sequence present in <seq>)
		H    Hard clip on the read (clipped sequence NOT present in <seq>)
		P    Padding (silent deletion from the padded reference sequence)
		 */
		boolean isAllM=false;
		String original;
		char[] alignedSeq;
		List<CigarElement> cigarElems;
		CigarOperator currentOp;//curent operation corresponding to the currentIndCigarElem (M,D,I,S,X,=)
		int basesLeft;
		int totalNewBases=0;
		int alignedBases=0;////index of new aligned bases in alignedSeq 
		int originalBases=0;//index of original bases in original sequence

		public CigarCounter(Cigar cigar,String readSeq){
			original =readSeq;
			cigarElems=cigar.getCigarElements();
			//check if the sequence needs to be altered (an all match cigar is going to have the same sequence before and after this class treatment)
			if(cigarElems.size()==1 && cigarElems.get(0).getOperator().equals(CigarOperator.M)){
				isAllM=true;
				System.out.println( "cigarElems: ALL Ms");
			}else{
				System.out.print( "cigarElems:");
				getTotalNewBases( cigar);
				alignedSeq=new char[totalNewBases];
				for (int e=0;e<cigarElems.size();e++){
					currentOp=cigarElems.get(e).getOperator();
					basesLeft=cigarElems.get(e).getLength();
					System.out.print( currentOp+" "+basesLeft+ " ");

					while(basesLeft>0){
						align();					
					}
				}
				//printout
				System.out.println( );
				System.out.println("alignedSeq:" );
				for (int i=0;i<alignedSeq.length;i++){
					System.out.print(alignedSeq[i] );
				}
				System.out.println( );
			}

		}

		private void getTotalNewBases(Cigar cigar) {
			for (int e=0;e<cigarElems.size();e++){
				if (cigar.getCigarElement(e).getOperator().equals(CigarOperator.M) || cigar.getCigarElement(e).getOperator().equals(CigarOperator.D) || cigar.getCigarElement(e).getOperator().equals(CigarOperator.EQ)  || cigar.getCigarElement(e).getOperator().equals(CigarOperator.I)){
					totalNewBases+=cigarElems.get(e).getLength();	
				}			
			}
		}

		private void align(){

			switch(currentOp){
			case M://Match					
				alignedSeq[alignedBases++]=original.charAt(originalBases++);
				break;
			case EQ://equal match
				alignedSeq[alignedBases++]=original.charAt(originalBases++);
				break;
			case D://deletion
				alignedSeq[alignedBases++]='-';
				break;
			case I://insertion: don't know if this is correct yet
				alignedSeq[alignedBases++]='*';
				break;
			case S://soft clipped (Ignore this base)
				originalBases++;
				break;
			case H://hard clipped (Ignore this base)
				originalBases++;
				break;
			case P://padding (empty insert) : don't know if this is correct yet
				alignedSeq[alignedBases++]='*';
				break;
			default://'X': sequence doesn't appear on alignment (don't map)->if position concerns the variation, signal to ignore the variation(returns -1?)
				;
				break;
			}
			basesLeft--;

		}

		private String getSubseq(int beg,int end){
			String subSeq="";
			if (!isAllM){
				for (int i=beg;i<end;i++){
					if(i<alignedSeq.length){
						subSeq+=alignedSeq[i];
					}else {
						subSeq="#";
						break;
					}
				}
			}else subSeq= original.substring(beg,end);
			return subSeq;
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
