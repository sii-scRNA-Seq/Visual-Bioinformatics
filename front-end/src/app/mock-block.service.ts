import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';

import { Block, BlockId } from './block.interface';
import { BlockServiceInterface } from './block.service.interface';

@Injectable({
  providedIn: 'root'
})
export class MockBlockService implements BlockServiceInterface {
  private readonly blocksOnCanvas$: BehaviorSubject<Block[]> = new BehaviorSubject<Block[]> ([]);
  readonly blocksOnCanvas: Observable<Block[]> = this.blocksOnCanvas$.asObservable();

  private readonly executingBlocks$: BehaviorSubject<boolean> = new BehaviorSubject<boolean> (false);
  readonly executingBlocks: Observable<boolean> = this.executingBlocks$.asObservable();

  addBlock(id: BlockId): void {
    id;
  }

  removeBlock(id: BlockId): void {
    id;
  }

  executeBlocks(): void {
    this.executingBlocks$.next(true);
  }
}
