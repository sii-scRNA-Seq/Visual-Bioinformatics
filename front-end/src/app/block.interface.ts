export interface Block {
  blockId: BlockId
  blockUUID: string
  title: string
  possibleChildBlocks: BlockId[]
  parameters: Parameter[]
}

export type BlockId = 'loaddata' | 'basicfiltering' | 'qcplots' | 'qcfiltering' | 'variablegenes' | 'pca' | 'integration' | 'runumap';
export const BlockIdToTitleMap: Record<BlockId, string> = {
  loaddata: 'Load Data',
  basicfiltering: 'Basic Filtering',
  qcplots: 'Quality Control Plots',
  qcfiltering: 'Quality Control Filtering',
  variablegenes: 'Identify Highly Variable Genes',
  pca: 'Principal Component Analysis',
  integration: 'Integration',
  runumap: 'Run UMAP'
};

export interface Parameter {
  type: 'InputParameter' | 'SelectParameter'
  key: string
  text: string
  value: number | string
  options?: Option[]
}

export interface Option {
  key: string
  text: string
}
