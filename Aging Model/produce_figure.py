from process import *
from plot import *


def main():
    ##Dataset 9 ######
    intersect_of_all = get_index_intersect(
    [
    'Oligo',
 'PC',
 'Astro',
 'OPC',
 'Inhibitory_Neur', 'Excitatory_Neur'], dataset_id = '9')

    list_of_predictions, y_test, ret = predict([
    'Oligo',
 'PC',
 'Astro',
 'OPC',
 'Inhibitory_Neur', 'Excitatory_Neur'],
    gene_num = 500, dataset_id='9', fixed_index = intersect_of_all)
    plot_table(ret, save_prefix=f'9')
    list_of_predictions, y_test, ret = predict([
    'Oligo',
 'PC',
 'Astro',
 'OPC',
 'Inhibitory_Neur', 'Excitatory_Neur'],
    gene_num = 2000, dataset_id='9', fixed_index = intersect_of_all)
    plot_table(ret, save_prefix=f'9')

    list_of_predictions, y_test, ret = predict( [
    'Oligo',
 'PC',
 'Astro',
 'OPC',
 'Inhibitory_Neur', 'Excitatory_Neur'],
    gene_num = 2000, dataset_id='9', fixed_index = None)
    plot_table(ret, save_prefix=f'9')
#### Dataset 24 #####

    intersect_of_all = get_index_intersect(
    [
    'L2.3.IT',
 'Oligo',
 'Astro',
 'L4.IT',
 'OPC', 'Chandelier__Pvalb', 'Sst__Sst.Chodl', 'L5.IT', 'L6.IT'], dataset_id = '24')
    list_of_predictions, y_test, ret = predict([
    'L2.3.IT',
 'Oligo',
 'Astro',
 'L4.IT',
 'OPC', 'Chandelier__Pvalb', 'Sst__Sst.Chodl', 'L5.IT', 'L6.IT'],
    gene_num = 2000, dataset_id='24', fixed_index = intersect_of_all)
    plot_table(ret, save_prefix=f'24')
    list_of_predictions, y_test, ret = predict([
    'L2.3.IT',
 'Oligo',
 'Astro',
 'L4.IT',
 'OPC', 'Chandelier__Pvalb', 'Sst__Sst.Chodl', 'L5.IT', 'L6.IT'],
    gene_num = 2000, dataset_id='24', fixed_index = None)
    plot_table(ret, save_prefix=f'24')

    ##### Dataset 21 #####
    intersect_of_all = get_index_intersect(
    [
    'Chandelier.Pvalb',
 'Astro',
 'L5_IT',
 'Vip',
 'L2.3_IT', 'Sst.Sst_Chodl', 'Oligo' , 'L4_IT', 'L6_IT'], dataset_id = '21')

    list_of_predictions, y_test, ret = predict([
    'Chandelier.Pvalb',
    'Astro',
    'L5_IT',
    'Vip',
    'L2.3_IT', 'Sst.Sst_Chodl', 'Oligo' , 'L4_IT', 'L6_IT'],
    gene_num = 500, dataset_id='21', fixed_index = intersect_of_all)

    plot_table(ret, save_prefix=f'21')

    list_of_predictions, y_test, ret = predict([
    'Chandelier.Pvalb',
    'Astro',
    'L5_IT',
    'Vip',
    'L2.3_IT', 'Sst.Sst_Chodl', 'Oligo' , 'L4_IT', 'L6_IT'],
    gene_num = 500, dataset_id='21', fixed_index = None)
    plot_table(ret, save_prefix=f'21')

if __name__ == '__main__':
    main()