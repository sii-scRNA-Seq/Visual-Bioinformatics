import { Observable } from 'rxjs';

import { Block } from './block.interface';
import { Output } from './output';

export interface OutputServiceInterface {
  outputs: Observable<Output[]>;

  executeBlocks(block: Block[]): void;
}
