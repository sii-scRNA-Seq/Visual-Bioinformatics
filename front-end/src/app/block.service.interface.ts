import { Observable } from 'rxjs';

import { Block, BlockId } from './block.interface';

export interface BlockServiceInterface {

  blocksOnCanvas: Observable<Block[]>;

  addBlock(id: BlockId): void;
  removeBlock(id: BlockId): void;
  executeBlocks(): void;
}
