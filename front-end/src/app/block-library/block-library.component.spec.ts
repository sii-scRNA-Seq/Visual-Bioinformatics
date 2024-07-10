import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';
import { MatTooltipModule } from '@angular/material/tooltip';

import { BlockId } from '../block.interface';
import { BlockLibraryComponent } from './block-library.component';
import { BlockService } from '../block.service';
import { MockBlockService } from '../mock-block.service';
import { MockOutputService } from '../mock-output.service';
import { OutputService } from '../output.service';

describe('BlockLibraryComponent', () => {
  let component: BlockLibraryComponent;
  let fixture: ComponentFixture<BlockLibraryComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [BlockLibraryComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
        MatTooltipModule
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: OutputService, useClass: MockOutputService },
      ],
    });
    fixture = TestBed.createComponent(BlockLibraryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Adding Blocks From Library', () => {
    it ('should call blockService.addBlock with corresponding block ID when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const blockIDs: BlockId[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      blockIDs.forEach(blockID => {
        const button = fixture.debugElement.query(By.css('#'.concat(blockID)));
        button.triggerEventHandler('click', {});
        fixture.detectChanges();
        expect(blockService.addBlock).toHaveBeenCalledWith(blockID);
      });
      expect(blockService.addBlock).toHaveBeenCalledTimes(blockIDs.length);
    });

    it ('should be available when blocks are not being executed', () => {
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('#loaddata')).nativeElement.disabled).toEqual(false);
    });

    it ('should become disabled while blocks are being executed', () => {
      component.executingBlocks = true;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const blockIDs: string[] = ['loaddata', 'basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      blockIDs.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(true);
      });
    });

    it ('should become available once blocks have stopped being executed', () => {
      component.executingBlocks = true;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('#loaddata')).nativeElement.disabled).toEqual(false);
    });

    it ('should only allow Load Data block to be added when canvas is empty', () => {
      component.blockList = [];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const enabledBlocks: string[] = ['loaddata'];
      const disabledBlocks: string[] = ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes', 'pca', 'integration', 'runumap'];
      enabledBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(false);
      });
      disabledBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(true);
      });
    });

    it ('should only allow possible child blocks of the last block on the canvas to be added when canvas is not empty', () => {
      component.blockList = [{
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'],
        parameters: [],
      }];
      component.executingBlocks = false;
      component.updateDisabledBlocks();
      fixture.detectChanges();
      const enabledBlocks: string[] = ['basicfiltering', 'qcplots', 'qcfiltering', 'variablegenes'];
      const disabledBlocks: string[] = ['loaddata', 'pca', 'integration', 'runumap'];
      enabledBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(false);
      });
      disabledBlocks.forEach(blockID => {
        expect(fixture.debugElement.query(By.css('#'.concat(blockID))).nativeElement.disabled).toEqual(true);
      });
    });
  });
});
