import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.security.KeyStore.Entry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
	List<Integer> tempReadIndexes = Arrays.asList( 269142, 269145, 269146, 269151, 269152, 269157, 269158, 269159, 269160); //ERASE after use. For debug only
	VCFParser vm;
	int varExpMatIndexesSize=0;
	HashMap<Integer,String> varExpMatIndexes=new HashMap<Integer,String>();//keeps track of the (possible) variation expressions index in the 
	//connectivity matrix (matIndex) <key> , and their id (VarId) <value>
	HashMap<String,VariationData> varExprIds=new HashMap<String,VariationData>();//keeps track of (possible)_variation_expresion_identifiers<key> and 
	//the corresponding (possible)_variation_expresion VariationData (which contains a matIndex field)
	List<VariationData> vcfVarLines = new ArrayList<VariationData>();//list of variations Positions where 'ref' different than 'alt'
	VariationsConnectivityMatrix vcm;//expressed variation connectivity matrix

	HashMap<Integer,ArrayList> readIndexVariationsIds=new HashMap <Integer,ArrayList> ();

	public VariationsManager(VCFParser vm){
		this.vm=vm;
	}




	public void registerPossibleVariationExpressions(VariationData vd){//for each variation position, two possible expressions are registered (Ref and alt)

		vd.setVarExpMatrixIndex(varExpMatIndexesSize++);
		varExpMatIndexes.put(vd.matrixIndex,vd.id);//stores reference	allele	
		varExprIds.put(vd.id, vd);
		//System.out.println(vd.id+"\t/"+vd.matrixIndex);

		VariationData newVD=vd.makeAlternativeVarExpressionFromRef(vd);	
		newVD.setVarExpMatrixIndex(varExpMatIndexesSize++);
		varExpMatIndexes.put(newVD.matrixIndex,newVD.id);//stores alternative allele	
		varExprIds.put(newVD.id, newVD);
		//if(newVD.id.equals("457C"))System.out.println("!!!!!!!!!!!!    457C REGISTRATION "+newVD.id+"\t/"+newVD.matrixIndex+" "+newVD.ref+" varExprIds.get:"+varExprIds.get(newVD.id).ref);
	}


	public void fillVariantMatrix() {

		//fill VariationsPos vector (vector containing all variations positions -kept vcf lines- ordered in a vector)
		int[] variationsPos=new int[vcfVarLines.size()];
		for (int v=0;v<vcfVarLines.size();v++){
			variationsPos[v]=vcfVarLines.get(v).pos;
		}	

		iterateSamFile(variationsPos);	
		//vcm.printMatrix(20);

	}

	private void iterateSamFile(int[] variationsPos) {
		SAMFileReader inputSam;
		int lastVarPos =variationsPos[variationsPos.length-1];System.out.print("lastVarPos:"+lastVarPos);

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
			CigarCounter cigCount = null;

			ArrayList<String> tempDetectedVariationsIds=new ArrayList<String>();
			int readIndex=0;//ERASE when FINISHED; For debug only

			while (iter.hasNext() /*&& ct<1500 */ ) {//iterates the sam file (for each read)
				SAMRecord rec = iter.next();//read line
				readName=rec.getReadName()+libraryName;
				alStart = rec.getAlignmentStart();
				alEnd = rec.getReadLength()+alStart;

				//identify  variation signatures of the read
				if (alStart<lastVarPos /*&& ct<1500 */){	//alStart<lastVarPos avoids any reads that doesn't starts before the last Variation position
					tempDetectedVariationsIds=new ArrayList<String>();

					//identify start var_position covered by read (first var):

					//chek that vpInd is correctly positioned	
					while(variationsPos[vpInd]<alStart)	vpInd++;//since bam is sorted, the index only increases
					currentVpInd=vpInd;//actualize current variant position index

					if(alEnd>variationsPos[currentVpInd]){
						try{
							readIndex=PloidyPhaser.readIndexes.get(readName);
						}catch (Exception e){
							readIndex=0;
							//System.err.println("Read Index error, readName:"+readName);
						}

						if (alStart>760 && alEnd<950)System.out.println("bam:"+libraryName+" ct:"+ct+" ReadIndex:"+PloidyPhaser.readIndexes.get(readName) +" st:"+alStart+" end:"+alEnd+" length:"+rec.getReadLength()+" cig:"+rec.getCigar()+" readSeq:"+rec.getReadString()+" variationsPos["+currentVpInd+"]:"+variationsPos[currentVpInd]);


						cigCount=new CigarCounter(rec.getCigar(),rec.getReadString(),readIndex);//new cigarCount for each read Sequence rec.getReadString()
						boolean endReached=false;
						//while the read spans positions with a variation
						while(!endReached   && alEnd>variationsPos[currentVpInd]  ){

							//get all vars covered by the read
							int i=0;//index over all posible variations for the current read
							VariationData curVar=vcfVarLines.get(currentVpInd+i);//currentVariation
							while(curVar.pos==variationsPos[currentVpInd] && !endReached){
								//curVar=vcfVarLines.get(currentVpInd+i);//currentVariation
								String subSeq=solveVariation(curVar, alStart,rec.getReadString(),rec.getReadLength(),cigCount);
								if (!subSeq.equals("#")){//if the subSequence has been solved
									tempDetectedVariationsIds.add(curVar.pos+subSeq);
									//System.out.println( "   CurrVar:"+curVar.outString()+" SUBSEQ "+subSeq);
								}
								i++;
								//System.out.print( "   i:"+i+" lastCurvar:"+curVar.id+" currentVpInd "+currentVpInd+ " variationsPos[currentVpInd]:"+variationsPos[currentVpInd]+" curVar.pos:"+curVar.pos+" variationsPos.length:"+variationsPos.length);
								if ((currentVpInd+i) <variationsPos.length){
									curVar=vcfVarLines.get(currentVpInd+i);
								}else endReached=true;;
							}
							currentVpInd=currentVpInd+i;
						}



					}
					//System.out.println(" end");
				}else break;//alStart > lastVarPos
				tempDetectedVariationsIds=checkDetectedVariations(tempDetectedVariationsIds,cigCount);
				if (alStart>760 && alEnd<950)System.out.println(" *read vars before check"+tempDetectedVariationsIds);
				//add this paired end variations to the readIndex variations (joint both paired end's variations in a single hashmap)
				if(readIndexVariationsIds.containsKey(readIndex) && readIndex!=0){//if the first paired end of the read has been registered. (readIndex=0 when there is a readIndex error)		
					for (int ndv=0;ndv<tempDetectedVariationsIds.size();ndv++){//add the variations from the second paired end of the read
						if (!readIndexVariationsIds.get(readIndex).contains(tempDetectedVariationsIds.get(ndv)))(readIndexVariationsIds.get(readIndex)).add(tempDetectedVariationsIds.get(ndv));//avoids repeats in overlaps paired ends 
					}
					//if (readIndex==437 || readIndex==438  || readIndex==459  || readIndex==463  )System.out.println(" read :"+readIndex+" final vars:"+readIndexVariationsIds.get(readIndex));

				}else if (readIndex!=0){//readIndex=0 when there is a readIndex error
					readIndexVariationsIds.put(readIndex, tempDetectedVariationsIds);
					if (alStart>760 && alEnd<950) System.out.println(" read new"+readIndex+" vars:"+tempDetectedVariationsIds+ " hasinsert?:"+cigCount.hasInsert);

				}


				ct++;
			}
			iter.close();
			inputSam.close();
		}

		//...and finally: fill the matrix
		Set<Map.Entry<Integer, ArrayList>> readIndexVariationsIdsSet = readIndexVariationsIds.entrySet();
		for(Map.Entry<Integer, ArrayList> ri : readIndexVariationsIdsSet){//for each readIndex (both pair ends together)
			int readInd =ri.getKey();
			ArrayList variationIDs=ri.getValue();


			for (int dvr=0;dvr<variationIDs.size();dvr++){//for each detected variation (row)
				for (int dvc=0;dvc<dvr;dvc++){//for each detected variation (column)
					int row=0;
					int col=0;
					try{
						row=varExprIds.get(variationIDs.get(dvr)).matrixIndex;
						col=varExprIds.get(variationIDs.get(dvc)).matrixIndex;

						vcm.addConection(row,col,readInd);
						//if (row==col)System.out.println	((dvr+dvc)+" dvr:"+dvr+" dvc:"+dvc+" readIndex:"+readIndex+" r==c ACCIDENT at row "+variationIDs.get(dvr)+" and column :"+variationIDs.get(dvc)+" r:"+row+" col:"+col+" variationIDs:"+variationIDs);

					}catch (Exception e){
						System.err.println(" Variation error at row "+variationIDs.get(dvr)+" and column :"+variationIDs.get(dvc)+" r:"+row+" col:"+col+" variationIDs:"+variationIDs);
					}
				}
			}
		}
	}


	private PairPosSignature getPairPosSignature(String id){

		StringBuilder p = new StringBuilder();
		StringBuilder s = new StringBuilder();
		for(int i = 0; i < id.length(); i++){
			char c = id.charAt(i);
			if(c > 47 && c < 58){
				p.append(c);
			}else s.append(c);
		}
		return new PairPosSignature(Integer.parseInt(p.toString()),s.toString());

	}

	private class PairPosSignature{
		int pos;
		String sig;

		private PairPosSignature(int p,String s){
			pos=p;
			sig=s;
		}
	}

	private ArrayList<String> checkDetectedVariations(ArrayList<String> tempDetectedVariationsIds, CigarCounter cigCount) {
		//System.out.println( " Temporary Detected VariationsIds:"+tempDetectedVariationsIds);
		ArrayList<String> correctedVariationsIds=new ArrayList<String>();
		String idTempCheckVar;
		int pos;
		String expSignature;
		int[]tempVarsToRemoveIndexes=new int[tempDetectedVariationsIds.size()];
		int tempVarsToRemove=0;

		for(int d=0;d<tempDetectedVariationsIds.size();d++){			
			idTempCheckVar=tempDetectedVariationsIds.get(d);			

			if(varExprIds.get(idTempCheckVar)==null){//we are interested in correcting variations if they have an error signature

				//System.out.print ("*******NEXT THAT NEEDS CHECK:  "+idTempCheckVar );
				PairPosSignature pps=getPairPosSignature(idTempCheckVar);//split pos from expressed variation part of the id
				pos=pps.pos;//position of var
				expSignature=pps.sig;//expressed allele

				//System.out.println (" **  id: "+idTempCheckVar+ " * :"+pos+"| |"+expSignature+": ***** ");


				if(expSignature.substring(0,1).equals("-")){//the pos is not covered by the alignment (an expressed deletion covers this position)-> remove it
					//System.out.println (idTempCheckVar+"  removed BECAUSE OF '-' " );
					tempVarsToRemoveIndexes[tempVarsToRemove++]=d;

				}else  if(expSignature.length()>1 && expSignature.substring(1,2).equals("-")){ //remove extra '-'

					//System.out.println (idTempCheckVar+"  (EXTRA -s)   Var Removed and replaced by "+pos+expSignature.substring(0,2) );

					tempVarsToRemoveIndexes[tempVarsToRemove++]=d;
					correctedVariationsIds.add(pos+expSignature.substring(0,2));

				} else { 
					//System.out.print ("   CHECKING COMBO-VARs FOR :" +idTempCheckVar);
					//TO DO : THIS PORTION SHOULD BE RECURSIVE FOR MORE RIGOUR (in long variations, there could also exist a sub position that needs to be corrected)

					//either is a combination of variations which signature can be corrected, or it is an error to dismiss

					//boolean isAVariationCombo=false;
					char[] correctedSig=expSignature.toCharArray();//the corrected version of the signature(initialize as the currently expressed)
					for (int esp=0;esp<expSignature.length();esp++){//for each sub-position of the expression

						//check existing variations in the corresponding positions
						int tp=pos+esp;// target position (starting at 0 )
						String qID=expSignature.substring(esp,esp+1);//candidate char in query signature
						String candSig=(tp+qID);//candidate signature
						int indOfCandVariation;//index of Candidate Variation in tempDetectedVariaton

						//if the candidate exists for this subposition
						if(tempDetectedVariationsIds.contains(candSig) ){
							//check the expression of the candiate subposition
							indOfCandVariation=tempDetectedVariationsIds.indexOf(candSig);//gets its index							
							//System.out.println ("   CANDIDATE COMBO-VAR candSig :" +candSig+/*": isAVariationComb?:"+isAVariationCombo+*/" at ind "+indOfCandVariation );

							if (varExprIds.containsKey(candSig)){
								correctedSig[esp]=varExprIds.get(candSig).ref.charAt(0);
								//System.out.print ("    CORRECTING "+idTempCheckVar+" at pos "+ esp+ " with  "+correctedSig[esp]); 
							}else{//the VariationCombo is not registered in the vcf file (probably a sequencng and/or alignment error)
								//System.out.println ("     BREAK LOOP for "+candSig+ "variation not registered in vcf file");

								break;
							}
						}

					}
					String correctedString=pos+new String(correctedSig);
					if (varExprIds.containsKey(correctedString)){//the correction has been succesful
						tempVarsToRemoveIndexes[tempVarsToRemove++]=d;//remove old
						correctedVariationsIds.add(correctedString);//replace by new
						//System.out.println("     CORRECTED Sig:"+correctedString);
					}else{//the VariationCombo is not registered in the vcf file
						tempVarsToRemoveIndexes[tempVarsToRemove++]=d;
						//System.out.println("###### CORRECTED Sig:"+correctedString+ " not found in varExprIds");
					}


				}

				//System.out.println (" ****************************" );


			} //else signature ok         
		}

		//System.out.println ("! tempDetectedVariationsIds:"+tempDetectedVariationsIds );
		//remove bad variations...
		for (int r=tempVarsToRemove-1;r>=0;r--){		
			tempDetectedVariationsIds.remove(tempVarsToRemoveIndexes[r]);
		}
		//System.out.println ("!! correctedVariationsIds:"+correctedVariationsIds );


		//...and add the corrected ones		
		if (correctedVariationsIds.size()>0){
			tempDetectedVariationsIds.addAll(correctedVariationsIds);
		}

		tempDetectedVariationsIds.sort(null);
		//System.out.println ("!!! tempDetectedVariationsIds:"+tempDetectedVariationsIds );
		cigCount=null;//to avoid using a cigCount from the wrong read when the corresponding loop is not entered in the parent method iterateSamFile()
		return tempDetectedVariationsIds;
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
			subSeq=cigCount.getSubseq(subSeqBeg, subSeqEnd, curVar);//get the aligned read subsequence corresponding to the pertinent positions
			if(curVar.pos==855){
				System.out.println( "   CurrVar:"+curVar.outString()+" GETSUBSEQ "+"("+subSeqBeg+"/"+subSeqEnd+"):"+subSeq);
				if (!cigCount.hasInsert && curVar.isInsert){

				}
			}

		}
		//System.out.println(" DETECTED "+"("+subSeqBeg+"/"+subSeqEnd+"):"+subSeq);

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



		public CigarCounter(Cigar cigar,String readSeq,int readIndex){
			original =readSeq;
			cigarElems=cigar.getCigarElements();
			//check if the sequence needs to be altered (an all match cigar is going to have the same sequence before and after this class treatment)
			if(cigarElems.size()==1 && cigarElems.get(0).getOperator().equals(CigarOperator.M)){//100M
				isAllM=true;
				//System.out.println( "     cigarElems: ALL Ms");
			}else{

				exploreCigar(cigar);
				alignedSeq=new char[totalNewBases];
				if (tempReadIndexes.contains(readIndex))System.out.print( "cigarElems:");


				for (int e=0;e<cigarElems.size();e++){					
					currentOp=cigarElems.get(e).getOperator();
					basesLeft=cigarElems.get(e).getLength();
					//if (tempReadIndexes.contains(readIndex))System.out.print( currentOp+" "+basesLeft+ " ");
					while(basesLeft>0 && isValidRead){
						align();					
					}

				}
				//printout
				/*
				if (isValidRead &&  (tempReadIndexes.contains(readIndex))){
					System.out.println( );
					System.out.print("alignedSeq:" );
					for (int i=0;i<alignedSeq.length;i++){
						System.out.print(alignedSeq[i] );
					}
					System.out.println( );
				}
				 */

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
				/*
				System.out.print("nbInsertRegions:"+nbInsertRegions +"  with lengths: ");
				for (int i=0;i<nbInsertRegions;i++){
					System.out.print(insertsRegionsLengths[i]+" pos:"+insertPositions[i]+", " );
				}
				System.out.println();
				 */
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

		private String getSubseq(int beg,int end, VariationData curVar){
			String subSeq="";
			
			boolean cigarHasInsertAtPos=false;
			//determine if the cigar has an insert at the position of the begining of this subseq
			if (!isAllM && hasInsert){
				for (int i=0;i<insertPositions.length;i++){
					if(insertPositions[i]==beg+1){
						cigarHasInsertAtPos=true;
						System.out.println("    --curVar:"+curVar.id+" --------   cigar HasInsertAtPos:"+(beg+1)+" original seq:"+original);
					}
				}
			}
			
			
			
			try{


				if (!isAllM ){//complex cigar
					int nextInsertPos=0;

					for (int i=beg;i<end;i++){
						if(!hasInsert){//cigar has no insert. 
							
							
							if(i<alignedSeq.length){
								subSeq+=alignedSeq[i];
							}else {
								subSeq="#";
								isValidRead=false;
								break;
							}
						}else{//if it has at least an insert 
							if(i<alignedSeq.length ){
								//System.out.print("nextInsertPos:"+nextInsertPos+ " insertPositions.length:"+insertPositions.length+" i:"+i);

								if(nextUpIndex!=alignedSeq.length)nextUpIndex=insertPositions[nextInsertPos];//updates nextUpIndex if the end of insertPositions has not been reached yet, otherwise, the value of alignedSeq.length has been set and is kept
								//System.out.println(" insertPositions[nextInsertPos]:"+insertPositions[nextInsertPos]+ " nextUpIndex:"+nextUpIndex+" addedBases:"+addedBases);

								if((i-addedBases)<nextUpIndex){//normal match of the sequence 
									subSeq+=alignedSeq[i+addedBases];
									//System.out.println(" -i:"+i+" add:"+alignedSeq[i+addedBases]);
								}else if ((i-addedBases)>=nextUpIndex){//insert match of the sequence
									//System.out.println("ENTER INSERT PART i:"+i+"addedBases++:"+addedBases);
									while(insertsRegionsLengths[nextInsertPos]>0){
										subSeq+=alignedSeq[(i-addedBases)];
										addedBases++;
										//System.out.println(" ++++addedBases++:"+addedBases);
										//System.out.println(" *i:"+(i-addedBases)+" add:"+alignedSeq[(i-addedBases)]);
										insertsRegionsLengths[nextInsertPos]--;
									} 

									if (insertPositions.length>(nextInsertPos+1)){
										nextInsertPos++;
									}else nextUpIndex=alignedSeq.length;//which runs the normal match of the sequence loop until the end of the read 

								}

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
		vcm= new VariationsConnectivityMatrix(varExpMatIndexes.size());

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

	public void colourMatrix(int minReadsThreshold) {//tries to cluster the readIndexes in k different groups. minReadsThreshold is the minimum number of reads to be cosidered when connecting two different variations
		System.out.println("colourMatrix ("+minReadsThreshold+")");

		//find Seed Nodes

		//count number of conexions per varation
		for (int r=0;r<vcm.matSize;r++){//for all variations (per row)
			int conectivityDegree=0;
			int c;
			System.out.print(varExpMatIndexes.get(r)+": ");

			for ( c=0;c<r;c++){//with the fixed row, move through the columns until r=c
				if(vcm.getReadList(r, c).size()>minReadsThreshold) {
					conectivityDegree+=vcm.getReadList(r, c).size();
					//System.out.print(" "+r+"/"+c+":"+vcm.getReadList(r, c).size());
				}
			}
			c=r;//now we fix the column and move through the rest of the rows 
			for (int rr=c;rr<vcm.matSize;rr++){
				if(vcm.getReadList(rr, c).size()>minReadsThreshold) {
					conectivityDegree+=vcm.getReadList(rr, c).size();
					//System.out.print(" "+rr+"/"+c+":"+vcm.getReadList(rr, c).size());
				}
			}
			System.out.println(conectivityDegree);

		}
	}
}
