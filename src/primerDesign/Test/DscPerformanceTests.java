package primerDesign.Test;

import java.text.NumberFormat;
import java.util.Stack;

import primerDesign.util.SimpleTimer;
import cern.colt.list.IntArrayList;

public class DscPerformanceTests {
	public static void main(String[] args){
		SimpleTimer timer = new SimpleTimer();
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		
		int elements = 100000;
		
		{
			System.out.println("Creating Integer stack, push & pop " + elements + " elements");
			Stack<Integer> stack = new Stack<Integer>();
			for(int i=0; i<elements; i++){
				stack.push(i);
			}
			System.out.println("Pushing took " + timer.getTimeString());
			for(int i=0; i<elements; i++){
				 stack.pop();
			}
			System.out.println("Popping took " + timer.getTimeString());
			System.gc();
			System.out.println("Heap: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
		
		{
			System.out.println("Creating IntArrayList, push & pop " + elements + " elements");
			IntArrayList list = new IntArrayList();
			int size = 0;
			for (int i = 0; i < elements; i++){
				list.add(i);
				size++;
			}
			System.out.println("Pushing took " + timer.getTimeString());
			for(int i=0; i<elements; i++){
				 list.get(i);
				 --size;
			}
			System.out.println("Popping took " + timer.getTimeString());
			System.gc();
			System.out.println("Heap: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
		
		{
			System.out.println("Creating int[], push & pop " + elements + " elements");
			int[] array = new int[0];
			int[] temp;
			for(int i=0; i<elements; i++){
				temp = new int[array.length+1];
				System.arraycopy(array, 0, temp, 0, array.length);
				temp[array.length] = i;
				array = temp;
			}
			System.out.println("Pushing took " + timer.getTimeString());
			int tmp;
			for(int i=0; i<array.length; i++){
				tmp = array[i];
			}
			System.out.println("Popping took " + timer.getTimeString());
			System.gc();
			System.out.println("Heap: " + format.format(runtime.totalMemory()-runtime.freeMemory()));
		}
	}
}
