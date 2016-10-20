import java.util.HashMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class VariationsManager {
	VCFManager vm;
	int varMatIndexesSize=0;
	HashMap<VariationData,Integer> varMatIndexes=new HashMap<VariationData,Integer>();//keeps track of vcf positions<key> that contains a 'ref' 
															//different than an 'alt', and their Index order <value>,in order to build the matrix,
	HashMap<String,VariationData> varDatIds=new HashMap<String,VariationData>();//keeps track of vcf identifiers<key> and the 
																				//corresponding VariationData 
	int[][]varMat;
	
	
	public VariationsManager(VCFManager vm){
		this.vm=vm;
	}
	
	
	public void registerVariation(VariationData vd){
		
		vd.setMatrixIndex(++varMatIndexesSize);
		varMatIndexes.put(vd,vd.matrixIndex);//stores reference	allele	
		varDatIds.put(vd.id, vd);
		System.out.print(vd.id+"/"+vd.matrixIndex);
		VariationData newVD=vd.makeAlterativeVarDataFromRef(vd);	
		newVD.setMatrixIndex(varMatIndexesSize);
		varMatIndexes.put(newVD,newVD.matrixIndex);//stores alternative allele	
		varDatIds.put(newVD.id, newVD);
		System.out.println("   "+newVD.id+"/"+newVD.matrixIndex);
		System.out.print("I THINK THESE TWO HASHMAP ARE OT CORRECT FOR WHAT I NEED... CHANGE TO GET TE INDEX, THE VARDAT AND THE ID CORRECTLY");
		
	}
	
	
	public void fillVariantMatrix() {
		
		SAMFileReader inputSam;
		for (int bi=0;bi<PloidyPhaser.bamPaths.size();bi++){
			inputSam = new SAMFileReader(PloidyPhaser.bamPaths.get(bi));





			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter = inputSam.iterator();
			//PrintWriter writer = new PrintWriter(Ploest.outputFile + "//" + Ploest.projectName+ "//"+Ploest.projectName+"SamParsed.txt", "UTF-8");
			//String line = "";
			int ct = 0;
			//System.out.println("Analyzing ");
			String refName = "";// for debugging a bad line in the .sam file
			int alStart;
			int alEnd;
			String libraryName="."+PloidyPhaser.bamPaths.get(bi).getName().substring(0,PloidyPhaser.bamPaths.get(bi).getName().lastIndexOf('.'));
			String readName="";
			
			while (iter.hasNext() && ct<10 ) {//iterates the sam file
				SAMRecord rec = iter.next();
				readName=rec.getReadName()+libraryName;
				refName = rec.getReferenceName();
				alStart = rec.getAlignmentStart();
				alEnd = rec.getReadLength()+alStart;
				if (alStart<10000 && ct<10 ){
					System.out.println("ct:"+ct+" read:"+readName+" ReadIndex:"+PloidyPhaser.readIndexes.get(readName) +" refName:"+refName+" st:"+alStart+" end:"+alEnd+" length:"+rec.getReadLength()+" cig:"+rec.getCigar()+" s:"+rec.getReadString());
					ct++;
				}
				//line = (i + " " + rec.getReadName() + " " + refName + " " + alStart + " " + rec.getAlignmentEnd() + " "
				//			+ rec.getReadLength() + " " + rec.getMappingQuality() + " ;");
				//writer.println(line);
			}
			iter.close();
			inputSam.close();
		}
	}
	
	
	public  void constructVariantMatrix() {
		varMat=new int[varMatIndexes.size()][varMatIndexes.size()];
		
	}
}
