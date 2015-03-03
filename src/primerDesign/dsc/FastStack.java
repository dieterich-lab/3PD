package primerDesign.dsc;

import weka.core.FastVector;

/**
 * This class implements a fast (un-synchronized!) integer stack.
 * 
 * @author Sebastian Fršhler
 *
 */
public class FastStack {
	FastVector stack;
	
	public FastStack(){
		stack = new FastVector();
	}
	
	public void push(int element){
		this.stack.addElement(element);
	}
	
	public int pop(){
		int result = (Integer) this.stack.lastElement();
		this.stack.removeElementAt(this.stack.size()-1);
		return result;
	}
	
	public int peek(){
		return (Integer) this.stack.lastElement();
	}
	
	public boolean empty(){
		return this.stack.size() == 0;
	}
	
	public String toString(){
		StringBuffer buffer = new StringBuffer();
		for(int i=0; i<stack.size(); i++){
			buffer.append(stack.elementAt(i));
			if(i<stack.size()-1) buffer.append("\t");
		}
		return buffer.toString();
	}
}
