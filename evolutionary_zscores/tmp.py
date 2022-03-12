THRESHOLD = 50

z_locations_s1 = pd.read_csv('zscores/NC_004718.3.csv', index_col=0)
z_locations_r = pd.read_csv('zscores/NC_045512.2.csv', index_col=0)
z_locations_s1['Sequence_id'] = 'NC_004718.3'
z_locations_r['Sequence_id'] = 'NC_045512.2'
z_locations_s1['Strain'] = 'Sars1'
z_locations_r['Strain'] = 'Reference'

z_locations_s1 = z_locations_s1[z_locations_s1['Z-Score'] >= THRESHOLD].reset_index(drop=True)
z_locations_r = z_locations_r[z_locations_r['Z-Score'] >= THRESHOLD].reset_index(drop=True)


base_offset = {
    'OL869974.1': 177,
    'NC_045512.2': 0,
}


z_locations_s1['Start'] += z_locations_s1.apply(
    lambda row: (
        aligned_sequences[row['Sequence_id']][base_offset[row['Sequence_id']] : row['Start']].count('-')
        + base_offset[row['Sequence_id']]
    ),
    axis=1,
)
z_locations_s1['End'] += z_locations_s1.apply(
    lambda row: (
        aligned_sequences[row['Sequence_id']][base_offset[row['Sequence_id']] : row['End']].count('-')
        + base_offset[row['Sequence_id']]
    ),
    axis=1,
)
