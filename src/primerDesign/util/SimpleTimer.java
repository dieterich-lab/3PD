package primerDesign.util;

import java.text.NumberFormat;

/**
 * Implements a simple timer for benchmarking purposes.
 * 
 * @author Sebastian Fršhler
 *
 */
public class SimpleTimer {

	private long time;
	private long totalTime;
	private NumberFormat format;
	
	public SimpleTimer(){
		this.startTimer();
		this.format = NumberFormat.getInstance();
	}
	
	public void startTimer(){
		time = System.currentTimeMillis();
		totalTime = System.currentTimeMillis();
	}
	
	public void resetTimer(){
		this.startTimer();
	}
	
	public long getTime(){
		long currentTime = (System.currentTimeMillis() - time);
		time = System.currentTimeMillis();
		return currentTime;
	}
	
	public String getTimeString(){
		long currentTime = this.getTime();
		return format.format(currentTime) + "ms";
	}
	
	public long getTotalTime(){
		return (System.currentTimeMillis() - totalTime);
	}
	
	public String getTotalTimestring(){
		return format.format(getTotalTime()) + "ms";
	}
}
