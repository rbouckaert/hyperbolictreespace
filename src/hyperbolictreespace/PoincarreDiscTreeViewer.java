package hyperbolictreespace;



import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.apache.commons.math3.util.FastMath;

import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.ClusterTree;
import beast.util.Randomizer;


@Description("Poincarre Disc tree:\n"
		+ "Distances are in hyperbolic geometric between points on circle\n"
		+ "Sliders determine angle of points on cirlce.")
public class PoincarreDiscTreeViewer extends CubeTreeViewer {
	private static final long serialVersionUID = 1L;

	public PoincarreDiscTreeViewer() {	
	}
	
	public PoincarreDiscTreeViewer(String startValues) {
		N = 4;
		initialise(startValues);
		
		animationState = new int[N];
		//Arrays.fill(animationState, 1);
	}
	
	protected void initialise(String startValues) {
		setLayout(new BorderLayout());
		controlPanel = new ControlPanel(startValues);
		add(controlPanel, BorderLayout.WEST);
		treePanel = new TreePanel2();
		add(treePanel, BorderLayout.CENTER);
	}
	

	public class TreePanel2 extends TreePanel {
		private static final long serialVersionUID = 1L;

		@Override
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			if (nodes == null) {
				return;
			}
			
			((Graphics2D)g).setStroke(new BasicStroke(1f));
			int R = getWidth();
			g.setColor(new Color(255,200,200));
			g.drawOval(R/4, R/4, R/2, R/2);
			g.setColor(Color.red);
			g.drawOval(0, 0, R, R);
			double h0 = 1.0;
			for (int i = 0; i < N; i++) {
				h0 = Math.min(h0, sliderR[i].getValue()/100.0);
			}
			System.out.println("h0 = " + h0);
			
			for (int i = 0; i < N; i++) {
				double angle1 = 2 * Math.PI * slider[i].getValue() / 360.0;
				double r = sliderR[i].getValue() / 100.0;
				double x1 = r * FastMath.sin(angle1);
				double y1 = r * FastMath.cos(angle1);
				g.fillOval((int)(R/2 + x1* R/2)-3, (int)(R/2 + y1*R/2)-3, 7, 7);
			}
			
			double [] angle = new double[nodes.length];
			for (int i = 0; i < N; i++) {
				angle[i] = slider[i].getValue();
			}
			for (int i = N; i < nodes.length; i++) {
				angle[i] = (angle[nodes[i].getLeft().getNr()] + angle[nodes[i].getRight().getNr()])/ 2.0; 
			}
			for (int i = 0; i < N; i++) {
				double angle1 = 2 * Math.PI * slider[i].getValue() / 360.0;
				double r = sliderR[i].getValue() / 100.0;
				double x1 = r * FastMath.sin(angle1);
				double y1 = r * FastMath.cos(angle1);

				Node parent = nodes[i].getParent();
				double angle2 = 2 * Math.PI * angle[parent.getNr()] / 360.0;
				r = 1.0-nodes[i].getParent().getHeight() / root.getHeight();
				r *= h0;
				double x2 = r * FastMath.sin(angle2);
				double y2 = r * FastMath.cos(angle2);

				g.drawLine((int)(R/2 + x1* R/2), (int)(R/2 + y1*R/2), (int)(R/2 + x2* R/2), (int)(R/2 + y2*R/2));
				
				g.drawString("" + (char) (65+i), (int)(R/2 + x1* R/2), (int)(R/2 + y1*R/2));
			}
			for (int i = N; i < nodes.length-1; i++) {
				double angle1 = 2 * Math.PI * angle[i] / 360.0;
				double r = 1.0-nodes[i].getHeight() / root.getHeight();
				r *= h0;
				double x1 = r * FastMath.sin(angle1);
				double y1 = r * FastMath.cos(angle1);

				Node parent = nodes[i].getParent();
				double angle2 = 2 * Math.PI * angle[parent.getNr()] / 360.0;
				r = 1.0-nodes[i].getParent().getHeight() / root.getHeight();
				r *= h0;
				double x2 = r * FastMath.sin(angle2);
				double y2 = r * FastMath.cos(angle2);

				g.drawLine((int)(R/2 + x1* R/2), (int)(R/2 + y1*R/2), (int)(R/2 + x2* R/2), (int)(R/2 + y2*R/2));
			}
			
		}
	}
	
	JSlider [] slider;
	JSlider [] sliderR;
	boolean useHyperBolic = true;

	public class ControlPanel extends JPanel {
		private static final long serialVersionUID = 1L;
		
		public ControlPanel(String startValues) {
			setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			
			slider = new JSlider[N];
			sliderR =  new JSlider[N];
			int [] start = new int[N];
			if (startValues != null) {
				String [] strs = startValues.split(",");
				for (int i = 0; i < N; i++) {
					start[i] = Integer.parseInt(strs[i]);
				}
			} else {
				for (int i = 0; i < N; i++) {
					start[i] = Randomizer.nextInt(360);
				}
			}
			for (int i = 0; i < N; i++) {
				slider[i] = new JSlider(0, 360, start[i]);
				slider[i].addChangeListener(new ChangeListener() {				
					@Override
					public void stateChanged(ChangeEvent e) {
						updateTree();
					}
				});
				
				JPanel sliderPanel = new JPanel();
				sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.X_AXIS));
				sliderPanel.add(new JLabel("" + (char)(65+i)));
				sliderPanel.add(slider[i]);
				add(sliderPanel);
			}

			for (int i = 0; i < N; i++) {
				sliderR[i] = new JSlider(0, 100, 50);
				sliderR[i].addChangeListener(new ChangeListener() {				
					@Override
					public void stateChanged(ChangeEvent e) {
						updateTree();
					}
				});
				
				JPanel sliderPanel = new JPanel();
				sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.X_AXIS));
				sliderPanel.add(new JLabel("R" + (char)(65+i)));
				sliderPanel.add(sliderR[i]);
				add(sliderPanel);
			}
			JCheckBox useHyperbolicCheckbox = new JCheckBox("<html>Use hyperbolic distance<br>instead of Euclidian</html>");
			useHyperbolicCheckbox.setSelected(useHyperBolic);
			useHyperbolicCheckbox.addActionListener(e -> {
				useHyperBolic = ((JCheckBox)e.getSource()).isSelected();
				updateTree();
			});
			add(useHyperbolicCheckbox);
		}
		

		
		public void updateTree() {
			int [] values = new int [N];
			double [] radius = new double [N];
			for (int i = 0; i < N; i++) {
				values[i] = slider[i].getValue();
				radius[i] = sliderR[i].getValue() / 100.0;
			}
			System.out.println(Arrays.toString(values));
			updateTree(values, radius);
		}
		


		public void updateTree(int [] values, double [] radius) {
			double [][] dist = new double[N][N];
			for (int i = 0; i < N; i++) {
				for (int j = i+1; j < N; j++) {
					dist[i][j] = distance(values[i], values[j], radius[i], radius[j], useHyperBolic);
					dist[j][i] = dist[i][j];
				}
			}
			
			Tree tree = matrix2Tree(dist, "upgma");
			root = tree.getRoot();
			nodes = tree.getNodesAsArray();
			
			System.out.print(root.toNewick());
			treePanel.repaint();
			
		}



	}
	
	public static Tree matrix2Tree(double[][] matrix, String clustertype) {
		int N = matrix.length;
		Distance distance = new Distance() {			
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				return matrix[taxon1][taxon2];
			}
		};
		
		TaxonSet taxonset = new TaxonSet();
		List<Sequence> sequences = new ArrayList<>();
		for (int i = 0; i < N; i++) {
			taxonset.taxonsetInput.get().add(new Taxon("taxon" + i));
			sequences.add(new Sequence("taxon" + i, "?"));
		}
		taxonset.initAndValidate();
		Alignment data = new Alignment(sequences, "nucleotide");
		
//		Matrix2Tree matrixTree = new Matrix2Tree();
//		matrixTree.initByName( 
//				"taxonset", taxonset,
//				"distance", distance);
//		return matrixTree;

		
		ClusterTree clusterTree = new ClusterTree();
		clusterTree.initByName("clusterType", clustertype, 
				"taxa", data,
				"distance", distance);
		return clusterTree;
	}
	
	/**
	 * 
	 * @param anglei in 
	 * @param anglej
	 * @param ri
	 * @param rj
	 * @param useHyperBolic
	 * @return
	 */
	public static double distance(double anglei, double anglej, double ri, double rj, boolean useHyperBolic) {
		double angle1 = 2 * Math.PI * anglei / 360.0;
		double x1 = ri * FastMath.sin(angle1);
		double y1 = ri * FastMath.cos(angle1);
		
		double angle2 = 2 * Math.PI * anglej / 360.0;
		double x2 = rj * FastMath.sin(angle2);
		double y2 = rj * FastMath.cos(angle2);
		double d = distance(angle1, angle2, ri, rj, x1, y1, x2, y2, useHyperBolic);
		return d;
	}
	
	/**
	 * 
	 * @param angle1
	 * @param angle2
	 * @param ri
	 * @param rj
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param useHyperBolic
	 * @return
	 */
	public static double distance(double angle1, double angle2,
			double ri, double rj,
			double x1, double y1,
			double x2, double y2,
			boolean useHyperBolic) {
		if (useHyperBolic) {
			double x = 2 * FastMath.sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1 - y2)) / ((1-ri) * (1-rj));
			double d = FastMath.acosh(1 + x);
			return d;
		}
		return FastMath.sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1 - y2));
	}
	
	
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.setSize(1024, 768);
		PoincarreDiscTreeViewer viewer = new PoincarreDiscTreeViewer(args.length > 0 ? args[0] : null);
		frame.add(viewer);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

}
