import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

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

		//fill Variations vector
		//System.out.println(" nbUniqueVars:"+nbUniqueVars+ " variations.size():"+variations.size());
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

			//System.out.println("Analyzing ");
			String refName = "";// for debugging a bad line in the .sam file
			int alStart;
			int alEnd;
			String libraryName="."+PloidyPhaser.bamPaths.get(bi).getName().substring(0,PloidyPhaser.bamPaths.get(bi).getName().lastIndexOf('.'));
			String readName="";
			int currentVpInd=0;//index of the current variation

			while (iter.hasNext() && ct<500 ) {//iterates the sam file (for each read)
				SAMRecord rec = iter.next();
				readName=rec.getReadName()+libraryName;
				refName = rec.getReferenceName();
				alStart = rec.getAlignmentStart();
				alEnd = rec.getReadLength()+alStart;
				//identify  variation signature of the read
				if (alStart<10000 && ct<500 ){					
					//identify start var_position covered by read (first var):
					
					//chek that vpInd is correctly positioned	
					while(variationsPos[vpInd]<alStart)	vpInd++;//since bam is sorted, the index only increases
					currentVpInd=vpInd;//actualize current variant position index
					if(alEnd>variationsPos[currentVpInd]){
						System.out.print("ct:"+ct+" ReadIndex:"+PloidyPhaser.readIndexes.get(readName) +" st:"+alStart+" end:"+alEnd+" length:"+rec.getReadLength()+" cig:"+rec.getCigar()+" s:"+rec.getReadString());
						System.out.println( "vpInd:"+vpInd+" currentVpInd:"+currentVpInd+" variationsPos_currentVpInd["+currentVpInd+"]:"+variationsPos[currentVpInd]);

					}
					
					//while the read spans positions with variation
					while(alEnd>variationsPos[currentVpInd]){
						//get all vars with same position
						int i=0;
						while(variations.get(currentVpInd+i).pos==variationsPos[currentVpInd]){
							System.out.println( " current variation:"+variations.get(currentVpInd+i).outString());
							//TODO for each variation check the ref or the alt expression and store the result in the matrix
							i++;

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
