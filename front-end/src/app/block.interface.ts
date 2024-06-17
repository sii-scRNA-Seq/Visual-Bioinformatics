export interface Block {
  blockId: BlockId
  title: string
  possibleChildBlocks: BlockId[]
  parameters: Parameter[]
}

export type BlockId = 'loaddata' | 'basicfiltering' | 'qcplots' | 'qcfiltering' | 'variablegenes' | 'pca' | 'runumap';

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
