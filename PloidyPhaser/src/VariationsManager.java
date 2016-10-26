import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
	HashMap<Integer,String> varMatIndexes=new HashMap<Integer,String>();//keeps track of the variation index in the 
	 							//connectivity matrix (matIndex) <key> , and their id (VarId) <value>
	HashMap<String,VariationData> varDatIds=new HashMap<String,VariationData>();//keeps track of variation_identifiers<key> and 
																				//the corresponding VariationData (which contains a matIndex field)
	List<VariationData> vcfVarLines = new ArrayList<VariationData>();//list of variations where 'ref' different than 'alt'
	int[][]varMat;


	public VariationsManager(VCFParser vm){
		this.vm=vm;
	}


	public void registerVariations(VariationData vd){

		vd.setMatrixIndex(varMatIndexesSize++);
		varMatIndexes.put(vd.matrixIndex,vd.id);//stores reference	allele	
		varDatIds.put(vd.id, vd);
		//System.out.println(vd.id+"\t/"+vd.matrixIndex);
		
		VariationData newVD=vd.makeAlterativeVarDataFromRef(vd);	
		newVD.setMatrixIndex(varMatIndexesSize++);
		varMatIndexes.put(newVD.matrixIndex,newVD.id);//stores alternative allele	
		varDatIds.put(newVD.id, newVD);
		//System.out.println(newVD.id+"\t/"+newVD.matrixIndex);
	}


	public void fillVariantMatrix() {

		//fill VariationsPos vector (vector containing all variations positions -kept vcf lines- ordered in a vector)
		int[] variationsPos=new int[vcfVarLines.size()];
		for (int v=0;v<vcfVarLines.size();v++){
			variationsPos[v]=vcfVarLines.get(v).pos;
		}	
/*
		//prints out varMatIndexes or varDatIds
		Iterator it = varMatIndexes.entrySet().iterator();
	    while (it.hasNext()) {
	    	Map.Entry pair = (Map.Entry)it.next();
	        System.out.println((pair.getKey()) + " = " + pair.getValue()+ " // "+varDatIds.get(pair.getValue()).matrixIndex);
	        //it.remove(); // avoids a ConcurrentModificationException
	    }
*/    
	    
		iterateSamFile(variationsPos);	
		
		
	}

	private void iterateSamFile(int[] variationsPos) {
		SAMFileReader inputSam;
		for (int bi=0;bi<PloidyPhaser.bamPaths.size();bi++){//for each bam file
			inputSam = new SAMFileReader(PloidyPhaser.bamPaths.get(bi));
			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter = inputSam.iterator();
			int ct = 0;

			int alStart;
			int alEnd;
			String libraryName="."+PloidyPhaser.bamPaths.get(bi).getName().substring(0,PloidyPhaser.bamPaths.get(bi).getName().lastIndexOf('.'));
			String readName="";
			int vpInd=0;//index pointing to the variationPos vector where the read maps its first variation
			int currentVpInd=0;//index of the current variation
			CigarCounter cigCount;
			List<String> tempDetectedVariationsIds=null;
			while (iter.hasNext() && ct<1500 ) {//iterates the sam file (for each read)
				SAMRecord rec = iter.next();//read line
				readName=rec.getReadName()+libraryName;
				alStart = rec.getAlignmentStart();
				alEnd = rec.getReadLength()+alStart;

				//identify  variation signatures of the read
				if (alStart<10000 && ct<1500 ){	
					tempDetectedVariationsIds=new ArrayList<String>();
					
					
					//identify start var_position covered by read (first var):

					//chek that vpInd is correctly positioned	
					while(variationsPos[vpInd]<alStart)	vpInd++;//since bam is sorted, the index only increases
					currentVpInd=vpInd;//actualize current variant position index
					if(alEnd>variationsPos[currentVpInd]){
						System.out.println("bam:"+libraryName+" ct:"+ct+" ReadIndex:"+PloidyPhaser.readIndexes.get(readName) +" st:"+alStart+" end:"+alEnd+" length:"+rec.getReadLength()+" cig:"+rec.getCigar()+" readSeq:"+rec.getReadString()+" variationsPos["+currentVpInd+"]:"+variationsPos[currentVpInd]);
						cigCount=new CigarCounter(rec.getCigar(),rec.getReadString());//new cigarCount for each read Sequence rec.getReadString()

						//while the read spans positions with a variation
						while(alEnd>variationsPos[currentVpInd]  ){

							//get all vars covered by the read
							int i=0;//index over all posible variations for the current read
							VariationData curVar=vcfVarLines.get(currentVpInd+i);//currentVariation
							while(curVar.pos==variationsPos[currentVpInd]){
								String subSeq=solveVariation(curVar, alStart,rec.getReadString(),rec.getReadLength(),cigCount);
								if (!subSeq.equals("#")){//if the subSequence has been solved
									tempDetectedVariationsIds.add(curVar.pos+subSeq);
									//System.out.println( "   CurrVar:"+curVar.outString()+" SUBSEQ "+subSeq);
									}
								
								//tempDetectedVariationsIds.add(solveVariation(curVar, alStart,rec.getReadString(),rec.getReadLength(),cigCount));
								i++;
								curVar=vcfVarLines.get(currentVpInd+i);
							}
							currentVpInd=currentVpInd+i;
						}
					}
				}
				System.out.println( "   tempDetectedVariationsIds:"+tempDetectedVariationsIds);
				ct++;
			}
			
			
			
			iter.close();
			inputSam.close();
		}
		
	}


	//for each variation, this function compares the subsequence of the candidate read to the possible allele variations
	//and outputs the represented one (ref=1 alt=2 none=0)
	public String solveVariation(VariationData curVar,int alStart, String readSeq, int readLength,CigarCounter cigCount){
		String subSeq="#";
		int subSeqBeg=curVar.pos-alStart;//begin index on the subsequence from readSeq.
		int subSeqEnd;//end index on the subsequence from readSeq.
		if (curVar.ref.length()>= curVar.alt.length()){
			subSeqEnd=(curVar.pos-alStart)+curVar.ref.length();//end index for SNP or deletion (Ref length >= alt)
		}else{
			subSeqEnd=(curVar.pos-alStart)+curVar.alt.length();//end index for insert (alt length > ref)
		}

		
		if(subSeqEnd<=readLength  && cigCount.isValidRead){	
			subSeq=cigCount.getSubseq(subSeqBeg, subSeqEnd);//get the aligned read subsequence corresponding to the pertinent positions
			//System.out.println( "   CurrVar:"+curVar.outString()+" GETSUBSEQ "+"("+subSeqBeg+"/"+subSeqEnd+"):"+subSeq);

		}
		return subSeq;


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
		boolean isValidRead=true;
		int[]insertsRegionsLengths;
		int[]insertPositions;
		boolean hasInsert=false;
		int addedBases=0;//keeps track of inserted bases to a getSubseq() request
		int nextUpIndex=0;//keeps track of upper bound index of sequence to getSubseq() request



		public CigarCounter(Cigar cigar,String readSeq){
			original =readSeq;
			cigarElems=cigar.getCigarElements();
			//check if the sequence needs to be altered (an all match cigar is going to have the same sequence before and after this class treatment)
			if(cigarElems.size()==1 && cigarElems.get(0).getOperator().equals(CigarOperator.M)){
				isAllM=true;
				System.out.println( "cigarElems: ALL Ms");
			}else{

				exploreCigar(cigar);
				alignedSeq=new char[totalNewBases];
				System.out.print( "cigarElems:");
				for (int e=0;e<cigarElems.size();e++){					
					currentOp=cigarElems.get(e).getOperator();
					basesLeft=cigarElems.get(e).getLength();
					System.out.print( currentOp+" "+basesLeft+ " ");
					while(basesLeft>0 && isValidRead){
						align();					
					}

				}
				//printout
				if (isValidRead){
					System.out.println( );
					System.out.print("alignedSeq:" );
					for (int i=0;i<alignedSeq.length;i++){
						System.out.print(alignedSeq[i] );
					}
					System.out.println( );
				}

			}

		}

		private void exploreCigar(Cigar cigar) {
			int nbInsertRegions=0;
			for (int e=0;e<cigarElems.size();e++){
				if (cigar.getCigarElement(e).getOperator().equals(CigarOperator.M) || cigar.getCigarElement(e).getOperator().equals(CigarOperator.D) || cigar.getCigarElement(e).getOperator().equals(CigarOperator.EQ)  || cigar.getCigarElement(e).getOperator().equals(CigarOperator.I)){
					totalNewBases+=cigarElems.get(e).getLength();	
				}	
				if(cigar.getCigarElement(e).getOperator().equals(CigarOperator.I)){	//computes number of Insert regions
					hasInsert=true;
					nbInsertRegions++;
				}				
			}
			if(hasInsert){//fills the vector inserts with each insert region length
				insertsRegionsLengths=new int[nbInsertRegions];
				insertPositions=new int[nbInsertRegions];
				int nextInsRegion=0;
				int normalBases=0;
				for (int e=0;e<cigarElems.size();e++){					
					if(cigar.getCigarElement(e).getOperator().equals(CigarOperator.I)){	
						insertsRegionsLengths[nextInsRegion]=cigar.getCigarElement(e).getLength();
						insertPositions[nextInsRegion++]=normalBases;
						normalBases=0;
					}else if (cigar.getCigarElement(e).getOperator().equals(CigarOperator.M) || cigar.getCigarElement(e).getOperator().equals(CigarOperator.D) || cigar.getCigarElement(e).getOperator().equals(CigarOperator.EQ)  ){
						normalBases+=cigar.getCigarElement(e).getLength();
					}
				}
				System.out.print("nbInsertRegions:"+nbInsertRegions +"  with lengths: ");
				for (int i=0;i<nbInsertRegions;i++){
					System.out.print(insertsRegionsLengths[i]+" pos:"+insertPositions[i]+", " );
				}
				System.out.println();
			}			
		}

		private void align(){
			try{
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
					hasInsert=true;
					alignedSeq[alignedBases++]=original.charAt(originalBases++);
					break;
				case S://soft clipped (Ignore this base)
					originalBases++;
					break;
				case H://hard clipped (Ignore this base)
					originalBases++;
					break;
				case P://padding (empty insert) : don't know if this is correct yet
					isValidRead=false;
					throw new ReadContainsException("Padding cigar operator ("+currentOp+") found in read "+original+"\nRead rejected");//alignedSeq[alignedBases++]='*';//throw new exception
					//break;
				default://'X': sequence doesn't appear on alignment (don't map)->if position concerns the variation, signal to ignore the variation(returns -1?)
					isValidRead=false;
					throw new ReadContainsException("unknown cigar operator ("+currentOp+") found in read "+original+"\nRead rejected");
				}
				basesLeft--;
			}catch (ReadContainsException e){
				System.out.println();
				System.out.println(e.getMessage());
			}
		}

		private String getSubseq(int beg,int end){
			String subSeq="";
			try{
			
			
			if (!isAllM ){
				int nextInsertPos=0;
				
				for (int i=beg;i<end;i++){
					if(!hasInsert){//cigar has no insert
						if(i<alignedSeq.length){
							subSeq+=alignedSeq[i];
						}else {
							subSeq="#";
							isValidRead=false;
							break;
						}
					}else{//if it has at least an insert
						if(i<alignedSeq.length){
							//System.out.print("nextInsertPos:"+nextInsertPos+ " insertPositions.length:"+insertPositions.length+" i:"+i);
							
						    if(nextUpIndex!=alignedSeq.length)nextUpIndex=insertPositions[nextInsertPos];//updates nextUpIndex if the end of insertPositions has not been reached yet, otherwise, the value of alignedSeq.length has been set and is kept
							//System.out.println(" insertPositions[nextInsertPos]:"+insertPositions[nextInsertPos]+ " nextUpIndex:"+nextUpIndex+" addedBases:"+addedBases);
							
							if((i-addedBases)<nextUpIndex){//normal match of the sequence 
								subSeq+=alignedSeq[i+addedBases];
								//System.out.println(" -i:"+i+" add:"+alignedSeq[i+addedBases]);
							}else if ((i-addedBases)==nextUpIndex){//insert match of the sequence
								while(insertsRegionsLengths[nextInsertPos]>0){
									subSeq+=alignedSeq[(i-addedBases)];
									addedBases++;
									//System.out.println(" *i:"+(i-addedBases)+" add:"+alignedSeq[(i-addedBases)]);
									insertsRegionsLengths[nextInsertPos]--;
								}
								if (insertPositions.length>(nextInsertPos+1)){
									nextInsertPos++;
								}else nextUpIndex=alignedSeq.length;//which runs the normal match of the sequence loop until the end of the read 

							}
							/*
							if (insertPositions.length>(nextInsertPos+1)){
								nextInsertPos++;
							}else nextUpIndex=alignedSeq.length;//which runs the normal match of the sequence loop until the end of the read 
*/
						}else {//index falls out of read
							subSeq="#";
							isValidRead=false;
							break;
						}
					}
				}
			}else subSeq= original.substring(beg,end);
			}catch (Exception e){
				System.err.print("Error in VariationsManager getSubseq("+ beg+", "+end+")" +" read sequence:"+original+" cigar:");
				for (int ce=0;ce<cigarElems.size();ce++){					
					currentOp=cigarElems.get(ce).getOperator();
					basesLeft=cigarElems.get(ce).getLength();
					System.err.print( currentOp+" "+basesLeft+ " ");
				}
				System.err.println();
			}
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
			for (int i = 0; i < (vcfVarLines.size()-1); i++) {
				System.out.println(vcfVarLines.get(i).outString());
			}
			System.out.print(vcfVarLines.get(vcfVarLines.size()-1).outString());//last line without endline
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


	class ReadContainsException extends Exception
	{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		//Parameterless Constructor
		public ReadContainsException() {}

		//Constructor that accepts a message
		public ReadContainsException(String message)
		{
			super(message);
		}
	}
}