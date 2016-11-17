
public class VariationData {
	String id="";//pos + either(reference or alternative allele)+ '-' before or after indicates an inseert or a deletion
   //The id expresses, not only  the vcf variation as such, but also which of the two alternative it expresses
	//id=pos+expSignature+(opt -)
	
	int pos;//position of the variation
	String expSignature;//expressed variation part of the signature
	String name="";//chromosome name
	
	String ref="";//the ref in the corresponding vcf line
	String alt=".";//the alt in the corresponding vcf line
	
	String sample="0/0";
	boolean isInsert=false;
	boolean isDeletion=false;
	boolean isBlank=false;
	int matrixIndex;


	public VariationData(){
		
	}
	public VariationData(String n,int p, String r,String a,String s ){
		name=n;
		pos=p;
		ref=r;
		alt=a;
		sample=s;
		id=Integer.toString(pos)+ref;//set regular id with ref
		expSignature=ref;
		if (ref.length()>alt.length())isDeletion=true;
		if (alt.length()>ref.length()){//is an insert
			isInsert=true;
			id+="-";//set id with insert identifier '-'
			expSignature+="-";
		}
	}
	
	public void printOut(){
		System.out.println(this.outString());
	}
	
	public String outString(){
		return /*name+"\t"+*/pos+"\t"+ref+"\t"+alt+"\t"/*+sample*/;		
	}
	
	public VariationData makeAlternativeVarExpressionFromRef(VariationData vd) {
		VariationData nvd=new VariationData(vd.name,vd.pos,vd.ref,vd.alt,vd.sample);
		nvd.id=Integer.toString(pos)+alt;//set regular id with alt
		nvd.expSignature=alt;
		if (ref.length()>alt.length()){//is a deletion
			isDeletion=true;
			nvd.id+="-";//set id with insert identifier '-'
			nvd.expSignature+="-";
		}
		return nvd;
	}
	
	public VariationData makeBlankVarExpressionFromRef(VariationData vd) {
		VariationData nvd=new VariationData(vd.name,vd.pos,vd.ref,vd.alt,vd.sample);
		nvd.isBlank=true;
		nvd.id=Integer.toString(pos)+"-";//set regular id with '-'
		nvd.expSignature="-";
		return nvd;
	}
	
	public void setVarExpMatrixIndex(int i){
		matrixIndex=i;
	}
	
	
}
