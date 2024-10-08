import pandas as pd
from kernel.Knormalized import knormalized
from kernel.subkernel import generative_sub_kernel
from kernel.kernel_alignment import alignment
from kernel.kernel_combination import drug_matrix_combination, atc_matrix_combination
from shortest_path.ATC_kernel import Similarity, Layered
from Network_consistency_projection.consistency_projection import NSP
from Top_Similar.WKNKN import WKNKN
import warnings
from data_split.cross_validation import Cross_validation
warnings.filterwarnings("ignore")


class Options(object):
    def __init__(self, drug_atc_path, atc_target_protein_path, atc_side_effects_path, atc_fingerprint_path, level, omega):
        drug_fingerprint_path, drug_side_effects_path, drug_target_protein_path, drug_interaction_path = self.file_path()
        self.drug_fingerprint = pd.read_csv(drug_fingerprint_path, index_col=0)
        self.drug_side_effects = pd.read_csv(drug_side_effects_path, index_col=0)
        self.drug_target_protein = pd.read_csv(drug_target_protein_path, index_col=0)
        self.drug_interaction = pd.read_csv(drug_interaction_path, index_col=0)
        self.drug_atc = pd.read_csv(drug_atc_path, index_col=0)
        self.atc_target_protein = pd.read_csv(atc_target_protein_path, index_col=0)
        self.atc_side_effects = pd.read_csv(atc_side_effects_path, index_col=0)
        self.atc_fingerprint = pd.read_csv(atc_fingerprint_path, index_col=0)
        self.level = level
        self.omega = omega


    def file_path(self):
        drug_fingerprint_path = 'data/drug_fingerprint/fingerprint.csv'
        drug_side_effects_path = './data/drug_side_effects/side_effects.csv'
        drug_target_protein_path = './data/drug_target_protein/uniprot.csv'
        drug_interaction_path = './data/drug_interaction/interaction_kernel.csv'
        return drug_fingerprint_path, drug_side_effects_path, drug_target_protein_path, drug_interaction_path

    def train(self, k):
        K_gip_fingerprint, K_corr_fingerprint, K_cos_fingerprint, K_jaccard_fingerprint, K_mi_fingerprint = generative_sub_kernel(self.drug_fingerprint.values)
        K_gip_effects, K_corr_effects, K_cos_effects, K_jaccard_effects, K_mi_effects = generative_sub_kernel(self.drug_side_effects.values)
        K_gip_protein, K_corr_protein, K_cos_protein, K_jaccard_protein, K_mi_protein = generative_sub_kernel(self.drug_target_protein.values)

        empty_predict_dataframe = pd.DataFrame()
        empty_actual_dataframe = pd.DataFrame()

        CV = Cross_validation(self.drug_atc)
        for positive_array, negative_array, train_data in CV.one2one_K_fold(k):

            K_gip_atc, K_corr_atc, K_cos_atc, K_jaccard_atc, K_mi_atc = generative_sub_kernel(train_data.values)
            K_ideal_drug = knormalized(train_data.values)
            drug_fingerprint_kernel = alignment(K_ideal_drug, K_gip_fingerprint, K_corr_fingerprint, K_cos_fingerprint, K_jaccard_fingerprint, K_mi_fingerprint)
            drug_side_effects_kernel = alignment(K_ideal_drug, K_gip_effects, K_corr_effects, K_cos_effects, K_jaccard_effects, K_mi_effects)
            drug_target_protein_kernel = alignment(K_ideal_drug, K_gip_protein, K_corr_protein, K_cos_protein, K_jaccard_protein, K_mi_protein)
            drug_atc_kernel = alignment(K_ideal_drug, K_gip_atc, K_corr_atc, K_cos_atc, K_jaccard_atc, K_mi_atc)
            drug_interaction_kernel = self.drug_interaction.values

            # 核融合
            drug_global_kernel = drug_matrix_combination(drug_fingerprint_kernel, drug_atc_kernel,
                                                         drug_interaction_kernel, drug_target_protein_kernel,
                                                         drug_side_effects_kernel)

            K_ideal_ATC = knormalized(train_data.values.T)

            K_gip_atc_effects, K_corr_atc_effects, K_cos_atc_effects, K_jaccard_atc_effects, K_mi_atc_effects = generative_sub_kernel(self.atc_side_effects.values)
            K_gip_atc_protein, K_corr_atc_protein, K_cos_atc_protein, K_jaccard_atc_protein, K_mi_atc_protein = generative_sub_kernel(self.atc_target_protein.values)
            K_gip_atc_fingerprint, K_corr_atc_fingerprint, K_cos_atc_fingerprint, K_jaccard_atc_fingerprint, K_mi_atc_fingerprint = generative_sub_kernel(self.atc_fingerprint.values)

            ATC_target_protein_kernel = alignment(K_ideal_ATC, K_gip_atc_protein, K_corr_atc_protein, K_cos_atc_protein, K_jaccard_atc_protein, K_mi_atc_protein)
            ATC_side_effects_kernel = alignment(K_ideal_ATC, K_gip_atc_effects, K_corr_atc_effects, K_cos_atc_effects, K_jaccard_atc_effects, K_mi_atc_effects)
            ATC_fingerprint_kernel = alignment(K_ideal_ATC, K_gip_atc_fingerprint, K_corr_atc_fingerprint, K_cos_atc_fingerprint, K_jaccard_atc_fingerprint, K_mi_atc_fingerprint)
            K_gip_ATC, K_corr_ATC, K_cos_ATC, K_jaccard_ATC, K_mi_ATC = generative_sub_kernel(train_data.values.T)
            ATC_atc_columns_kernel = alignment(K_ideal_ATC, K_gip_ATC, K_corr_ATC, K_cos_ATC, K_jaccard_ATC, K_mi_ATC)


            SPro = Similarity(train_data)
            ATC_atc_probabilistic_kernel = SPro.probabilistic_kernel(
                path=f'./shortest_path/new_{self.level}ATC_shortest_path_length_matrix.csv')

            SM = Layered(list(self.drug_atc.columns), level=self.level)
            ATC_atc_SM_kernel = SM.get_SM_kernel()

            atc_global_kernel = atc_matrix_combination(
                [ATC_atc_columns_kernel, ATC_atc_probabilistic_kernel, ATC_atc_SM_kernel, ATC_target_protein_kernel, ATC_side_effects_kernel, ATC_fingerprint_kernel])

            WK = WKNKN(drug_global_kernel, atc_global_kernel, train_data.values, self.omega)
            F_train_inference = WK.get_scores()


            nsp = NSP(drug_global_kernel, atc_global_kernel, F_train_inference)
            predict = nsp.network_NSP()

            positive_index = positive_array[:, 0]
            positive_column = positive_array[:, 1]
            positive_scores = pd.DataFrame(predict[positive_index, positive_column])

            negative_index = negative_array[:, 0]
            negative_column = negative_array[:, 1]
            negative_scores = pd.DataFrame(predict[negative_index, negative_column])

            empty_predict_dataframe = pd.concat([empty_predict_dataframe, positive_scores])
            empty_predict_dataframe = pd.concat([empty_predict_dataframe, negative_scores])

            actual_positive = pd.DataFrame(self.drug_atc.values[positive_index, positive_column])
            actual_negative = pd.DataFrame(self.drug_atc.values[negative_index, negative_column])

            empty_actual_dataframe = pd.concat([empty_actual_dataframe, actual_positive])
            empty_actual_dataframe = pd.concat([empty_actual_dataframe, actual_negative])


        eval_predict = empty_predict_dataframe
        eval_actual = empty_actual_dataframe
        eval_predict.to_csv(f'TEST_omega_{self.omega}_level_{self.level}_predict.csv', index=False)
        eval_actual.to_csv(f'TEST_omega_{self.omega}_level_{self.level}_actual.csv', index=False)


if __name__ == "__main__":
    drug_atc_path = 'data/drug_ATC/second_ATC.csv'
    atc_target_protein_path = 'data/ATC_target_protein/second_uniprot.csv'
    atc_side_effects_path = 'data/ATC_side_effects/second_side_effects.csv'
    atc_fingerprint_path = 'data/ATC_fingerprint/second_fingerprint.csv.csv'
    op = Options(drug_atc_path, atc_target_protein_path, atc_side_effects_path, atc_fingerprint_path, 2, 0.9)
    op.train(10)
