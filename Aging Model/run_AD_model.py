from process import *


def parse_argument():
    # argumant parser
    parser = argparse.ArgumentParser(description='Process data from a given csv file and return a dataframe')
    parser.add_argument('-g', '--gene_num', help='Number of genes', type=int, default=-1)
    parser.add_argument('-d', '--dataset_id', help='Dataset id', default='24')
    parser.add_argument('-p', '--partition', help = 'partition',type=int,default=1)
    parser.add_argument('-i', '--index', help = 'index of the partition',type=int, default=0)
    parser.add_argument('-s', '--stratify', help = 'Stratified sampling', type=int, default=0)
    parser.add_argument('-c', '--cross_val', help = 'Cross Validation Fold', type = int, default=0)
    parser.add_argument('-n', '--num_runs', help = 'Number of different seeds', type = int, default=1)
    parser.add_argument('-A', '--AD_Model', help= 'Train on control and predict on AD', type = int, default=1)
    parser.add_argument('--save_intermediate', help='save intermediate results', type = int, default=0)
    parser.add_argument('--shap', help='perform shap', type = int, default=0)
    parser.add_argument('--feature_importance', help='record feature importance', type = int, default=0)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_argument()
    celltype_names =  [
    'L2.3.IT',
    'Oligo',
    'Astro',
    'L4.IT',
    'OPC', 'Chandelier__Pvalb', 'Sst__Sst.Chodl' , 'L5.IT', 'L6.IT', 'Vip', 'Micro.PVM']



    # celltype_names = ['L2.3.IT']

    # celltype_names = ['L2.3.IT']
    # intersect_of_all = get_index_intersect(celltype_names, dataset_id = args.dataset_id)
    intersect_of_all = None
    celltype_names_subsetted = [celltype_names[i] for i in range(args.index, len(celltype_names), args.partition)]
   
    rets = predict_with_different_seeds(celltype_names_subsetted,
    gene_num = args.gene_num, dataset_id=args.dataset_id, 
    fixed_index = intersect_of_all, 
    AD_model = False, stratify = args.stratify, 
    grid_search_cv = args.cross_val, 
    num_splits = args.num_runs,
    AD_Model = args.AD_Model, 
    save_intermediate = args.save_intermediate, 
    shap = args.shap, 
    feature_importance = args.feature_importance)

    stratify = 'stratified' if args.stratify else 'random'
    rets.to_csv(f'results_{args.index}_outof_{args.partition}_gene_num_{args.gene_num}_dataset_{args.dataset_id}_{stratify}.csv')