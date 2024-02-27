export interface Block {
  blockId: BlockId
  title: string
  possibleChildBlocks: BlockId[]
  parameters: Parameter[]
}

export type BlockId = 'loaddata' | 'basicfiltering' | 'qcplots' | 'qcfiltering' | 'variablegenes' | 'pca';

export interface Parameter {
  key: string
  text: string
  value: number
}
  