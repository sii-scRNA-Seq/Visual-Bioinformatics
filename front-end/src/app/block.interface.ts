import { Observable } from "rxjs"

export interface Block {
    blockId: BlockId
    title: string
    possibleChildBlocks: BlockId[]
    parameters: Parameter[]
    onRun: (block: Block) => Observable<unknown>
}

export type BlockId = 'loaddata' | 'basicfiltering';

export interface Parameter {
    key: string
    text: string
    value: number
}
  