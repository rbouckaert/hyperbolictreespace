package hyperbolictreespace;




import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import javax.imageio.ImageIO;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.util.FastMath;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.app.util.TreeFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.matrixalgebra.IllegalDimension;

@Description("Convert posterior tree set to Poincarre disk using hyperbolic geometry.")
public class Trees2PoincarreDisk extends beast.core.Runnable {

	final public Input<TreeFile> treesInput = new Input<>("trees", "beast.trees on which this operation is performed",
			new TreeFile("[[none]]"));
	final public Input<TreeFile> mccInput = new Input<>("mcc", "summary tree used to guide the taxon order",
			new TreeFile("[[none]]"));
	final public Input<OutFile> outputInput = new Input<>("out", "output file. Fails if not specified", Validate.REQUIRED);
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
			"percentage of trees to used as burn-in (and will be ignored). Set >= 100 if only a cube is required", 10);

	final public Input<OutFile> pngInput = new Input<>("png", "output file for Poincarre disk visualisation. Ignored if not specified", new OutFile("[[none]]"));

	final static int imageSize = 1024;


	// maps taxon name to ordering
	private Map<String, Integer> taxonOrder;
	// taxa in order of taxonOrder
	private String[] taxa;
	private Color[] taxonColour;

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		TreeFileParser parser0 = new TreeFileParser(mccInput.get().getAbsolutePath(), 0);
		Tree [] mcctrees = parser0.parseFile();

		// sanity check
		if (mcctrees.length > 1) {
			throw new IllegalArgumentException("Expected only a single tree in tree file " + mccInput.get().getAbsolutePath());
		}
		// determine leaf order compatible with post order traversal of MCC tree
		taxonOrder = new HashMap<>();
		order(mcctrees[0].getRoot(), new int[1]);
		taxa = new String[taxonOrder.size()];
		for (String taxon : taxonOrder.keySet()) {
			taxa[taxonOrder.get(taxon)] = taxon;
		}
		taxonColour = new Color[taxa.length];
		for (int i = 0; i < taxa.length; i++) {
			taxonColour[i] = Color.getHSBColor((float)i/taxa.length, 0.9f, 0.5f);
		}
		for (int i = 0; i < taxa.length; i++) {
			System.err.println(taxa[i]);
		}
		
		// load tree set
		TreeFileParser parser = new TreeFileParser(treesInput.get().getAbsolutePath(), burnInPercentageInput.get());
		Tree [] trees = parser.parseFile();
		Log.warning(trees.length + " trees fetched");
		String [] taxa2 = trees[0].getTaxaNames();

		// sanity check
		if (!checkTaxaAreCompatibl(taxa, taxa2)) {
			throw new IllegalArgumentException("Taxon sets in MCC tree and tree set are different");
		}

		process(trees);

		Log.warning("Done");
	}

	private boolean checkTaxaAreCompatibl(String[] taxa, String[] taxa2) {
		if (taxa.length != taxa2.length) {
			// different number of taxa
			return false;
		}
		for (String taxon : taxa) {
			if (indexof(taxon, taxa2) < 0) {
				return false;
			}
		}
		return true;
	}

	private int indexof(String taxon, String[] taxa2) {
		for (int i = 0; i < taxa2.length; i++) {
			if (taxon.equals(taxa2[i])) {
				return i;
			}
		}
		return -1;
	}

	private void order(Node node, int [] i) {
		if (node.isLeaf()) {
			taxonOrder.put(node.getID(), i[0]++);
		} else {
			order(node.getLeft(), i);
			order(node.getRight(), i);
		}
	}



	private void process(Tree [] trees) throws IllegalDimension, IOException {
		BufferedImage bi;
		Graphics2D g;
		bi = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_RGB);
		g = (Graphics2D) bi.getGraphics();
		g.setPaintMode();
		g.setColor(Color.white);
		g.fillRect(0, 0, imageSize, imageSize);
		g.setColor(Color.blue);
		
		g.transform(AffineTransform.getScaleInstance(1, -1));
		g.transform(AffineTransform.getTranslateInstance(0, -imageSize) );
		
		g.drawLine(imageSize/2-5, imageSize/2, imageSize/2+5, imageSize/2);
		g.drawLine(imageSize/2, imageSize/2-5, imageSize/2, imageSize/2+5);
		int delta = imageSize*3/4 / trees[0].getLeafNodeCount(); 
		for (int i = 1; i <= taxa.length; i++) {
			g.drawLine(i * delta-5, imageSize/4, i * delta+5, imageSize/4);
			g.drawLine(i * delta, imageSize/4-5, i * delta, imageSize/4+5);			
		}

		
		normaliseTrees(trees);
		
		
		
		double scale = 1;
		int k = 0;
		for (Tree tree : trees) {
			double [][] points = tree2PoinCarreDisk(tree, taxonOrder);

			int i = 0;
			for (double [] point  : points) {
				g.setColor(taxonColour[i++]);
				point[0] *= scale;
				point[1] *= scale;
				int x = (int)(point[0] * imageSize/2 + imageSize/2);
				int y = (int)(point[1] * imageSize/2 + imageSize/2);
				g.drawOval(x, y, 1, 1);

				x = (int)(point[0] * imageSize/2 + i * delta);
				y = (int)(point[1] * imageSize/2 + imageSize/4);
				g.drawOval(x, y, 1, 1);
				
				
				if (i < points.length) {
					g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.025f));
				    int x2 = (int)(points[i][0] * scale * imageSize/2 + (i+1) * delta);
					int y2 = (int)(points[i][1] * scale * imageSize/2 + imageSize/4);
					g.drawLine(x, y, x2, y2);
					g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1.0f));
				}

				x = (int)(point[2] * scale * imageSize/2 + imageSize*2/3);
				y = (int)(point[3]/2 * imageSize/2 + imageSize/2);
				g.drawOval(x, y, 1, 1);

				x = (int)(point[3]/(2*Math.PI) * imageSize/2 + i * delta);
				y = (int)(point[2] * imageSize/3 + imageSize*6/8);
				g.drawOval(x, y, 1, 1);
				if (i < points.length) {
					g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.025f));
					int x2 = (int)(points[i][3] /(2*Math.PI) * imageSize/2 + (i+1) * delta);
					int y2 = (int)(points[i][2] * imageSize/3 + imageSize*6/8);
					g.drawLine(x, y, x2, y2);
					g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1.0f));
				}
			}
			
			k++;
		}

		ImageIO.write(bi, "png", pngInput.get());

	}



	private void normaliseTrees(Tree[] trees) {
		double h = trees[0].getRoot().getHeight();
		for (Tree tree : trees) {
			if (tree.getRoot().getHeight() > h) {
				h = tree.getRoot().getHeight();
			}
		}
		for (Tree tree : trees) {
			for (Node node : tree.getNodesAsArray()) {
				node.setHeight(node.getHeight() / h);
			}
		}		
	}

	private double[][] tree2PoinCarreDisk(Tree tree, Map<String, Integer> taxonOrder) {
		sort(tree.getRoot(), taxonOrder);
		// report(tree.getRoot());
		double [][] points= new double[taxonOrder.size()][4];
		double h = tree.getRoot().getHeight();
		double r = Math.tanh(h / 2);		
		traverse(tree.getRoot(), r, points, taxonOrder, new double[]{0.0});
		return points;
	}


	final static DecimalFormat f = new DecimalFormat("#.###");
	
	private void report(Node node, PrintStream out) {
		if (node.isLeaf()) {
			out.print(node.getID());
		} else {
			out.print("(");
			report(node.getLeft(), out);
			out.print(":");
			//out.print(f.format(node.getLeft().getLength()));
			out.print(",");
			report(node.getRight(), out);
			out.print(":");
			//out.print(f.format(node.getRight().getLength()));
			out.print(")");
		}
		if (node.isRoot()) {
			out.println();
		}		
	}

	
	/** recursively determines points of the tips of the tree 
	 * @param node current node in tree
	 * @param index keeps track of next order index to add
	 * @param internalNodes nodes for which heights are recorded
	 */
	private void traverse(Node node, double r, double[][] points, Map<String, Integer> taxonOrder,
			double [] angle) {
		if (node.isLeaf()) {
			int i = taxonOrder.get(node.getID());
			points[i][0] = r * Math.cos(angle[0]);
			points[i][1] = r * Math.sin(angle[0]);
			points[i][2] = r;
			points[i][3] = angle[0];
		} else {
			Node first = null, second = null;
			first = node.getChild(0);
			second = node.getChild(1);
			traverse(first, r, points, taxonOrder, angle);
			
			double distance = node.getHeight() * 2.0;
			angle[0] += angle(distance, r);
			//angle[0] = angle(distance, r);
			traverse(second, r, points, taxonOrder, angle);
		}		
	}
	
//	private void traverse(Node node, double r, double[][] points, Map<String, Integer> taxonOrder,
//			double [] angle, boolean [] visited) {
//		// System.out.println(node.toNewick(false) + " " + Arrays.toString(visited));
//		if (node.isLeaf()) {
//			if (Double.isInfinite(angle[0])) {
//				angle[0] = 0;
//				int i = taxonOrder.get(node.getID());
//				points[i][0] = r;
//				points[i][1] = 0;
//				points[i][2] = r;
//				points[i][3] = 0;
//				return;
//			}
//			Node g = node.getParent();
//			while (!g.isRoot() && visited[g.getNr()]) {
//				g = g.getParent();
//			}
//			visited[g.getNr()] = true;
//			double distance = g.getHeight() * 2.0;
//			
//			System.out.println(node.getID() + " " + distance/2);
//
//			// update angle so that newAngle = oldAngle + distance between this leaf and gap to previous leaf
//			// TODO: make sure this makes sense
//			angle[0] += angle(distance, r);
//			// angle[0] = angle(distance, r);
//
//			int i = taxonOrder.get(node.getID());
//			points[i][0] = r * Math.cos(angle[0]);
//			points[i][1] = r * Math.sin(angle[0]);
//			points[i][2] = r;
//			points[i][3] = angle[0];
//		} else {
//			traverse(node.getLeft(), r, points, taxonOrder, angle, visited);
//			traverse(node.getRight(), r, points, taxonOrder, angle, visited);
//		}
//	}

	private double angle(double distance, double r) {
		double x1 = 0;
		double y1 = r;
		
        UnivariateFunction f = new UnivariateFunction() {

			@Override
			public double value(double angle) {
				double x2 = r * FastMath.sin(angle);
				double y2 = r * FastMath.cos(angle);
				double x = 2 * FastMath.sqrt((x1-x2) * (x1-x2) + (y1-y2) * (y1 - y2)) / ((1-r) * (1-r));
				double d = FastMath.acosh(1 + x);
				double diff = Math.abs(distance - d);
				// System.err.println(angle + " " + diff);
				return diff;
			}
        };
		
        BrentOptimizer optimizer = new BrentOptimizer(1e-6, 1e-5);
        UnivariatePointValuePair p = optimizer.optimize(new MaxEval(200),
                new UnivariateObjectiveFunction(f),
                GoalType.MINIMIZE,
                new InitialGuess(new double [] {0.1}),
                new SearchInterval(0, Math.PI));
		double optimal = p.getPoint();
		return optimal;
	}

	
	private int sort(Node node, Map<String, Integer> taxonOrder) {
        if (node.isLeaf()) {
            return taxonOrder.get(node.getID());
        }

        int left = sort(node.getLeft(), taxonOrder);
        int right = sort(node.getRight(), taxonOrder);

        if (left < right) {
        	return left;
        }

        Node leftNode = node.getLeft();
        node.removeChild(leftNode);
        node.addChild(leftNode);
        return right;
	}

	
	
	public static void main(String[] args) throws Exception {
		new Application(new Trees2PoincarreDisk(), "Trees2PoincarreDisk", args);
	}


}
