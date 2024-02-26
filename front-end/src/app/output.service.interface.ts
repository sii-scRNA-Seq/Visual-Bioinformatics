import { Observable } from 'rxjs';

import { Block } from './block.interface';
import { Output } from './output';

export interface OutputServiceInterface {
  outputs: Observable<Output[]>;

  executeBlock(block: Block): Promise<void>;
}
