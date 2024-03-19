import { Observable } from 'rxjs';

import { Block, BlockId } from './block.interface';

export interface BlockServiceInterface {
  blocksOnCanvas: Observable<Block[]>;
  executingBlocks: Observable<boolean>;

  addBlock(id: BlockId): void;
  removeBlock(id: BlockId): void;
  executeBlocks(): void;
}
