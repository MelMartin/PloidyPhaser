
public class VariationData {
	String name;
	int pos;
	String ref;
	String alt;
	String sample;

	
	public VariationData(String n,int p, String r,String a,String s ){
		name=n;
		pos=p;
		ref=r;
		alt=a;
		sample=s;
	}
	
	public void printOut(){
		System.out.println(name+"\t"+pos+"\t"+ref+"\t:"+alt+"\t"+sample);
	}
}
